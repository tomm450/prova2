#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
a5foam = OpenFOAMReader(FileName='TO_BE_LOAD')
a5foam.SkipZeroTime = 1
a5foam.CaseType = 'Reconstructed Case'
a5foam.LabelSize = '32-bit'
a5foam.ScalarSize = '64-bit (DP)'
a5foam.Createcelltopointfiltereddata = 1
a5foam.Adddimensionalunitstoarraynames = 0
a5foam.MeshRegions = ['internalMesh']
a5foam.CellArrays = ['U', 'k', 'nut', 'omega', 'p']
a5foam.PointArrays = []
a5foam.LagrangianArrays = []
a5foam.Cachemesh = 1
a5foam.Decomposepolyhedra = 1
a5foam.ListtimestepsaccordingtocontrolDict = 0
a5foam.Lagrangianpositionswithoutextradata = 1
a5foam.Readzones = 0

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on a5foam
a5foam.MeshRegions = ['airfoil']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [887, 548]

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.LockDataRange = 0
pLUT.InterpretValuesAsCategories = 0
pLUT.ShowCategoricalColorsinDataRangeOnly = 0
pLUT.RescaleOnVisibilityChange = 0
pLUT.EnableOpacityMapping = 0
pLUT.RGBPoints = [-49100.94921875, 0.231373, 0.298039, 0.752941, -20495.27197265625, 0.865003, 0.865003, 0.865003, 8110.4052734375, 0.705882, 0.0156863, 0.14902]
pLUT.UseLogScale = 0
pLUT.ColorSpace = 'Diverging'
pLUT.UseBelowRangeColor = 0
pLUT.BelowRangeColor = [0.0, 0.0, 0.0]
pLUT.UseAboveRangeColor = 0
pLUT.AboveRangeColor = [1.0, 1.0, 1.0]
pLUT.NanColor = [1.0, 1.0, 0.0]
pLUT.Discretize = 1
pLUT.NumberOfTableValues = 256
pLUT.ScalarRangeInitialized = 1.0
pLUT.HSVWrap = 0
pLUT.VectorComponent = 0
pLUT.VectorMode = 'Magnitude'
pLUT.AllowDuplicateScalars = 1
pLUT.Annotations = []
pLUT.ActiveAnnotatedValues = []
pLUT.IndexedColors = []

# show data in view
a5foamDisplay = Show(a5foam, renderView1)
# trace defaults for the display properties.
a5foamDisplay.Representation = 'Surface'
a5foamDisplay.AmbientColor = [1.0, 1.0, 1.0]
a5foamDisplay.ColorArrayName = ['POINTS', 'p']
a5foamDisplay.DiffuseColor = [1.0, 1.0, 1.0]
a5foamDisplay.LookupTable = pLUT
a5foamDisplay.MapScalars = 1
a5foamDisplay.InterpolateScalarsBeforeMapping = 1
a5foamDisplay.Opacity = 1.0
a5foamDisplay.PointSize = 2.0
a5foamDisplay.LineWidth = 1.0
a5foamDisplay.Interpolation = 'Gouraud'
a5foamDisplay.Specular = 0.0
a5foamDisplay.SpecularColor = [1.0, 1.0, 1.0]
a5foamDisplay.SpecularPower = 100.0
a5foamDisplay.Ambient = 0.0
a5foamDisplay.Diffuse = 1.0
a5foamDisplay.EdgeColor = [0.0, 0.0, 0.5]
a5foamDisplay.BackfaceRepresentation = 'Follow Frontface'
a5foamDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
a5foamDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
a5foamDisplay.BackfaceOpacity = 1.0
a5foamDisplay.Position = [0.0, 0.0, 0.0]
a5foamDisplay.Scale = [1.0, 1.0, 1.0]
a5foamDisplay.Orientation = [0.0, 0.0, 0.0]
a5foamDisplay.Origin = [0.0, 0.0, 0.0]
a5foamDisplay.Pickable = 1
a5foamDisplay.Texture = None
a5foamDisplay.Triangulate = 0
a5foamDisplay.NonlinearSubdivisionLevel = 1
a5foamDisplay.UseDataPartitions = 0
a5foamDisplay.OSPRayUseScaleArray = 0
a5foamDisplay.OSPRayScaleArray = 'p'
a5foamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
a5foamDisplay.Orient = 0
a5foamDisplay.OrientationMode = 'Direction'
a5foamDisplay.SelectOrientationVectors = 'U'
a5foamDisplay.Scaling = 0
a5foamDisplay.ScaleMode = 'No Data Scaling Off'
a5foamDisplay.ScaleFactor = 0.12460310012102127
a5foamDisplay.SelectScaleArray = 'p'
a5foamDisplay.GlyphType = 'Arrow'
a5foamDisplay.UseGlyphTable = 0
a5foamDisplay.GlyphTableIndexArray = 'p'
a5foamDisplay.UseCompositeGlyphTable = 0
a5foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
a5foamDisplay.SelectionCellLabelBold = 0
a5foamDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
a5foamDisplay.SelectionCellLabelFontFamily = 'Arial'
a5foamDisplay.SelectionCellLabelFontSize = 18
a5foamDisplay.SelectionCellLabelItalic = 0
a5foamDisplay.SelectionCellLabelJustification = 'Left'
a5foamDisplay.SelectionCellLabelOpacity = 1.0
a5foamDisplay.SelectionCellLabelShadow = 0
a5foamDisplay.SelectionPointLabelBold = 0
a5foamDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
a5foamDisplay.SelectionPointLabelFontFamily = 'Arial'
a5foamDisplay.SelectionPointLabelFontSize = 18
a5foamDisplay.SelectionPointLabelItalic = 0
a5foamDisplay.SelectionPointLabelJustification = 'Left'
a5foamDisplay.SelectionPointLabelOpacity = 1.0
a5foamDisplay.SelectionPointLabelShadow = 0
a5foamDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
a5foamDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Arrow' selected for 'GlyphType'
a5foamDisplay.GlyphType.TipResolution = 6
a5foamDisplay.GlyphType.TipRadius = 0.1
a5foamDisplay.GlyphType.TipLength = 0.35
a5foamDisplay.GlyphType.ShaftResolution = 6
a5foamDisplay.GlyphType.ShaftRadius = 0.03
a5foamDisplay.GlyphType.Invert = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
a5foamDisplay.DataAxesGrid.XTitle = 'X Axis'
a5foamDisplay.DataAxesGrid.YTitle = 'Y Axis'
a5foamDisplay.DataAxesGrid.ZTitle = 'Z Axis'
a5foamDisplay.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
a5foamDisplay.DataAxesGrid.XTitleBold = 0
a5foamDisplay.DataAxesGrid.XTitleItalic = 0
a5foamDisplay.DataAxesGrid.XTitleFontSize = 12
a5foamDisplay.DataAxesGrid.XTitleShadow = 0
a5foamDisplay.DataAxesGrid.XTitleOpacity = 1.0
a5foamDisplay.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
a5foamDisplay.DataAxesGrid.YTitleBold = 0
a5foamDisplay.DataAxesGrid.YTitleItalic = 0
a5foamDisplay.DataAxesGrid.YTitleFontSize = 12
a5foamDisplay.DataAxesGrid.YTitleShadow = 0
a5foamDisplay.DataAxesGrid.YTitleOpacity = 1.0
a5foamDisplay.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
a5foamDisplay.DataAxesGrid.ZTitleBold = 0
a5foamDisplay.DataAxesGrid.ZTitleItalic = 0
a5foamDisplay.DataAxesGrid.ZTitleFontSize = 12
a5foamDisplay.DataAxesGrid.ZTitleShadow = 0
a5foamDisplay.DataAxesGrid.ZTitleOpacity = 1.0
a5foamDisplay.DataAxesGrid.FacesToRender = 63
a5foamDisplay.DataAxesGrid.CullBackface = 0
a5foamDisplay.DataAxesGrid.CullFrontface = 1
a5foamDisplay.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.ShowGrid = 0
a5foamDisplay.DataAxesGrid.ShowEdges = 1
a5foamDisplay.DataAxesGrid.ShowTicks = 1
a5foamDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
a5foamDisplay.DataAxesGrid.AxesToLabel = 63
a5foamDisplay.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
a5foamDisplay.DataAxesGrid.XLabelBold = 0
a5foamDisplay.DataAxesGrid.XLabelItalic = 0
a5foamDisplay.DataAxesGrid.XLabelFontSize = 12
a5foamDisplay.DataAxesGrid.XLabelShadow = 0
a5foamDisplay.DataAxesGrid.XLabelOpacity = 1.0
a5foamDisplay.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
a5foamDisplay.DataAxesGrid.YLabelBold = 0
a5foamDisplay.DataAxesGrid.YLabelItalic = 0
a5foamDisplay.DataAxesGrid.YLabelFontSize = 12
a5foamDisplay.DataAxesGrid.YLabelShadow = 0
a5foamDisplay.DataAxesGrid.YLabelOpacity = 1.0
a5foamDisplay.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
a5foamDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
a5foamDisplay.DataAxesGrid.ZLabelBold = 0
a5foamDisplay.DataAxesGrid.ZLabelItalic = 0
a5foamDisplay.DataAxesGrid.ZLabelFontSize = 12
a5foamDisplay.DataAxesGrid.ZLabelShadow = 0
a5foamDisplay.DataAxesGrid.ZLabelOpacity = 1.0
a5foamDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
a5foamDisplay.DataAxesGrid.XAxisPrecision = 2
a5foamDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
a5foamDisplay.DataAxesGrid.XAxisLabels = []
a5foamDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
a5foamDisplay.DataAxesGrid.YAxisPrecision = 2
a5foamDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
a5foamDisplay.DataAxesGrid.YAxisLabels = []
a5foamDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
a5foamDisplay.DataAxesGrid.ZAxisPrecision = 2
a5foamDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
a5foamDisplay.DataAxesGrid.ZAxisLabels = []

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
a5foamDisplay.PolarAxes.Visibility = 0
a5foamDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
a5foamDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
a5foamDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
a5foamDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
a5foamDisplay.PolarAxes.EnableCustomRange = 0
a5foamDisplay.PolarAxes.CustomRange = [0.0, 1.0]
a5foamDisplay.PolarAxes.PolarAxisVisibility = 1
a5foamDisplay.PolarAxes.RadialAxesVisibility = 1
a5foamDisplay.PolarAxes.DrawRadialGridlines = 1
a5foamDisplay.PolarAxes.PolarArcsVisibility = 1
a5foamDisplay.PolarAxes.DrawPolarArcsGridlines = 1
a5foamDisplay.PolarAxes.NumberOfRadialAxes = 0
a5foamDisplay.PolarAxes.AutoSubdividePolarAxis = 1
a5foamDisplay.PolarAxes.NumberOfPolarAxis = 0
a5foamDisplay.PolarAxes.MinimumRadius = 0.0
a5foamDisplay.PolarAxes.MinimumAngle = 0.0
a5foamDisplay.PolarAxes.MaximumAngle = 90.0
a5foamDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
a5foamDisplay.PolarAxes.Ratio = 1.0
a5foamDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.PolarAxisTitleVisibility = 1
a5foamDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
a5foamDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
a5foamDisplay.PolarAxes.PolarLabelVisibility = 1
a5foamDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
a5foamDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
a5foamDisplay.PolarAxes.RadialLabelVisibility = 1
a5foamDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
a5foamDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
a5foamDisplay.PolarAxes.RadialUnitsVisibility = 1
a5foamDisplay.PolarAxes.ScreenSize = 10.0
a5foamDisplay.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
a5foamDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
a5foamDisplay.PolarAxes.PolarAxisTitleBold = 0
a5foamDisplay.PolarAxes.PolarAxisTitleItalic = 0
a5foamDisplay.PolarAxes.PolarAxisTitleShadow = 0
a5foamDisplay.PolarAxes.PolarAxisTitleFontSize = 12
a5foamDisplay.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
a5foamDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
a5foamDisplay.PolarAxes.PolarAxisLabelBold = 0
a5foamDisplay.PolarAxes.PolarAxisLabelItalic = 0
a5foamDisplay.PolarAxes.PolarAxisLabelShadow = 0
a5foamDisplay.PolarAxes.PolarAxisLabelFontSize = 12
a5foamDisplay.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
a5foamDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
a5foamDisplay.PolarAxes.LastRadialAxisTextBold = 0
a5foamDisplay.PolarAxes.LastRadialAxisTextItalic = 0
a5foamDisplay.PolarAxes.LastRadialAxisTextShadow = 0
a5foamDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
a5foamDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
a5foamDisplay.PolarAxes.EnableDistanceLOD = 1
a5foamDisplay.PolarAxes.DistanceLODThreshold = 0.7
a5foamDisplay.PolarAxes.EnableViewAngleLOD = 1
a5foamDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
a5foamDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
a5foamDisplay.PolarAxes.PolarTicksVisibility = 1
a5foamDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
a5foamDisplay.PolarAxes.TickLocation = 'Both'
a5foamDisplay.PolarAxes.AxisTickVisibility = 1
a5foamDisplay.PolarAxes.AxisMinorTickVisibility = 0
a5foamDisplay.PolarAxes.ArcTickVisibility = 1
a5foamDisplay.PolarAxes.ArcMinorTickVisibility = 0
a5foamDisplay.PolarAxes.DeltaAngleMajor = 10.0
a5foamDisplay.PolarAxes.DeltaAngleMinor = 5.0
a5foamDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
a5foamDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
a5foamDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
a5foamDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
a5foamDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
a5foamDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
a5foamDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
a5foamDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
a5foamDisplay.PolarAxes.ArcMajorTickSize = 0.0
a5foamDisplay.PolarAxes.ArcTickRatioSize = 0.3
a5foamDisplay.PolarAxes.ArcMajorTickThickness = 1.0
a5foamDisplay.PolarAxes.ArcTickRatioThickness = 0.5
a5foamDisplay.PolarAxes.Use2DMode = 0
a5foamDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
a5foamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToNext()

animationScene1.GoToLast()

# create a new 'Slice'
slice1 = Slice(Input=a5foam)
slice1.SliceType = 'Plane'
slice1.Crinkleslice = 0
slice1.Triangulatetheslice = 1
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.37698449939489365, -0.02500000037252903, -0.03879009559750557]
slice1.SliceType.Normal = [1.0, 0.0, 0.0]
slice1.SliceType.Offset = 0.0

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, -0.01, 0.0]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, -0.01, 0.0]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.AmbientColor = [1.0, 1.0, 1.0]
slice1Display.ColorArrayName = ['POINTS', 'p']
slice1Display.DiffuseColor = [1.0, 1.0, 1.0]
slice1Display.LookupTable = pLUT
slice1Display.MapScalars = 1
slice1Display.InterpolateScalarsBeforeMapping = 1
slice1Display.Opacity = 1.0
slice1Display.PointSize = 2.0
slice1Display.LineWidth = 1.0
slice1Display.Interpolation = 'Gouraud'
slice1Display.Specular = 0.0
slice1Display.SpecularColor = [1.0, 1.0, 1.0]
slice1Display.SpecularPower = 100.0
slice1Display.Ambient = 0.0
slice1Display.Diffuse = 1.0
slice1Display.EdgeColor = [0.0, 0.0, 0.5]
slice1Display.BackfaceRepresentation = 'Follow Frontface'
slice1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
slice1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
slice1Display.BackfaceOpacity = 1.0
slice1Display.Position = [0.0, 0.0, 0.0]
slice1Display.Scale = [1.0, 1.0, 1.0]
slice1Display.Orientation = [0.0, 0.0, 0.0]
slice1Display.Origin = [0.0, 0.0, 0.0]
slice1Display.Pickable = 1
slice1Display.Texture = None
slice1Display.Triangulate = 0
slice1Display.NonlinearSubdivisionLevel = 1
slice1Display.UseDataPartitions = 0
slice1Display.OSPRayUseScaleArray = 0
slice1Display.OSPRayScaleArray = 'p'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.Orient = 0
slice1Display.OrientationMode = 'Direction'
slice1Display.SelectOrientationVectors = 'U'
slice1Display.Scaling = 0
slice1Display.ScaleMode = 'No Data Scaling Off'
slice1Display.ScaleFactor = 0.12460310012102127
slice1Display.SelectScaleArray = 'p'
slice1Display.GlyphType = 'Arrow'
slice1Display.UseGlyphTable = 0
slice1Display.GlyphTableIndexArray = 'p'
slice1Display.UseCompositeGlyphTable = 0
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelBold = 0
slice1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
slice1Display.SelectionCellLabelFontFamily = 'Arial'
slice1Display.SelectionCellLabelFontSize = 18
slice1Display.SelectionCellLabelItalic = 0
slice1Display.SelectionCellLabelJustification = 'Left'
slice1Display.SelectionCellLabelOpacity = 1.0
slice1Display.SelectionCellLabelShadow = 0
slice1Display.SelectionPointLabelBold = 0
slice1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
slice1Display.SelectionPointLabelFontFamily = 'Arial'
slice1Display.SelectionPointLabelFontSize = 18
slice1Display.SelectionPointLabelItalic = 0
slice1Display.SelectionPointLabelJustification = 'Left'
slice1Display.SelectionPointLabelOpacity = 1.0
slice1Display.SelectionPointLabelShadow = 0
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Arrow' selected for 'GlyphType'
slice1Display.GlyphType.TipResolution = 6
slice1Display.GlyphType.TipRadius = 0.1
slice1Display.GlyphType.TipLength = 0.35
slice1Display.GlyphType.ShaftResolution = 6
slice1Display.GlyphType.ShaftRadius = 0.03
slice1Display.GlyphType.Invert = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice1Display.DataAxesGrid.XTitle = 'X Axis'
slice1Display.DataAxesGrid.YTitle = 'Y Axis'
slice1Display.DataAxesGrid.ZTitle = 'Z Axis'
slice1Display.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.XTitleFontFamily = 'Arial'
slice1Display.DataAxesGrid.XTitleBold = 0
slice1Display.DataAxesGrid.XTitleItalic = 0
slice1Display.DataAxesGrid.XTitleFontSize = 12
slice1Display.DataAxesGrid.XTitleShadow = 0
slice1Display.DataAxesGrid.XTitleOpacity = 1.0
slice1Display.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.YTitleFontFamily = 'Arial'
slice1Display.DataAxesGrid.YTitleBold = 0
slice1Display.DataAxesGrid.YTitleItalic = 0
slice1Display.DataAxesGrid.YTitleFontSize = 12
slice1Display.DataAxesGrid.YTitleShadow = 0
slice1Display.DataAxesGrid.YTitleOpacity = 1.0
slice1Display.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
slice1Display.DataAxesGrid.ZTitleBold = 0
slice1Display.DataAxesGrid.ZTitleItalic = 0
slice1Display.DataAxesGrid.ZTitleFontSize = 12
slice1Display.DataAxesGrid.ZTitleShadow = 0
slice1Display.DataAxesGrid.ZTitleOpacity = 1.0
slice1Display.DataAxesGrid.FacesToRender = 63
slice1Display.DataAxesGrid.CullBackface = 0
slice1Display.DataAxesGrid.CullFrontface = 1
slice1Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.ShowGrid = 0
slice1Display.DataAxesGrid.ShowEdges = 1
slice1Display.DataAxesGrid.ShowTicks = 1
slice1Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
slice1Display.DataAxesGrid.AxesToLabel = 63
slice1Display.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.XLabelFontFamily = 'Arial'
slice1Display.DataAxesGrid.XLabelBold = 0
slice1Display.DataAxesGrid.XLabelItalic = 0
slice1Display.DataAxesGrid.XLabelFontSize = 12
slice1Display.DataAxesGrid.XLabelShadow = 0
slice1Display.DataAxesGrid.XLabelOpacity = 1.0
slice1Display.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.YLabelFontFamily = 'Arial'
slice1Display.DataAxesGrid.YLabelBold = 0
slice1Display.DataAxesGrid.YLabelItalic = 0
slice1Display.DataAxesGrid.YLabelFontSize = 12
slice1Display.DataAxesGrid.YLabelShadow = 0
slice1Display.DataAxesGrid.YLabelOpacity = 1.0
slice1Display.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
slice1Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
slice1Display.DataAxesGrid.ZLabelBold = 0
slice1Display.DataAxesGrid.ZLabelItalic = 0
slice1Display.DataAxesGrid.ZLabelFontSize = 12
slice1Display.DataAxesGrid.ZLabelShadow = 0
slice1Display.DataAxesGrid.ZLabelOpacity = 1.0
slice1Display.DataAxesGrid.XAxisNotation = 'Mixed'
slice1Display.DataAxesGrid.XAxisPrecision = 2
slice1Display.DataAxesGrid.XAxisUseCustomLabels = 0
slice1Display.DataAxesGrid.XAxisLabels = []
slice1Display.DataAxesGrid.YAxisNotation = 'Mixed'
slice1Display.DataAxesGrid.YAxisPrecision = 2
slice1Display.DataAxesGrid.YAxisUseCustomLabels = 0
slice1Display.DataAxesGrid.YAxisLabels = []
slice1Display.DataAxesGrid.ZAxisNotation = 'Mixed'
slice1Display.DataAxesGrid.ZAxisPrecision = 2
slice1Display.DataAxesGrid.ZAxisUseCustomLabels = 0
slice1Display.DataAxesGrid.ZAxisLabels = []

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.Visibility = 0
slice1Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
slice1Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
slice1Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
slice1Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
slice1Display.PolarAxes.EnableCustomRange = 0
slice1Display.PolarAxes.CustomRange = [0.0, 1.0]
slice1Display.PolarAxes.PolarAxisVisibility = 1
slice1Display.PolarAxes.RadialAxesVisibility = 1
slice1Display.PolarAxes.DrawRadialGridlines = 1
slice1Display.PolarAxes.PolarArcsVisibility = 1
slice1Display.PolarAxes.DrawPolarArcsGridlines = 1
slice1Display.PolarAxes.NumberOfRadialAxes = 0
slice1Display.PolarAxes.AutoSubdividePolarAxis = 1
slice1Display.PolarAxes.NumberOfPolarAxis = 0
slice1Display.PolarAxes.MinimumRadius = 0.0
slice1Display.PolarAxes.MinimumAngle = 0.0
slice1Display.PolarAxes.MaximumAngle = 90.0
slice1Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
slice1Display.PolarAxes.Ratio = 1.0
slice1Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.PolarAxisTitleVisibility = 1
slice1Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
slice1Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
slice1Display.PolarAxes.PolarLabelVisibility = 1
slice1Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
slice1Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
slice1Display.PolarAxes.RadialLabelVisibility = 1
slice1Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
slice1Display.PolarAxes.RadialLabelLocation = 'Bottom'
slice1Display.PolarAxes.RadialUnitsVisibility = 1
slice1Display.PolarAxes.ScreenSize = 10.0
slice1Display.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.PolarAxisTitleOpacity = 1.0
slice1Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
slice1Display.PolarAxes.PolarAxisTitleBold = 0
slice1Display.PolarAxes.PolarAxisTitleItalic = 0
slice1Display.PolarAxes.PolarAxisTitleShadow = 0
slice1Display.PolarAxes.PolarAxisTitleFontSize = 12
slice1Display.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.PolarAxisLabelOpacity = 1.0
slice1Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
slice1Display.PolarAxes.PolarAxisLabelBold = 0
slice1Display.PolarAxes.PolarAxisLabelItalic = 0
slice1Display.PolarAxes.PolarAxisLabelShadow = 0
slice1Display.PolarAxes.PolarAxisLabelFontSize = 12
slice1Display.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
slice1Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
slice1Display.PolarAxes.LastRadialAxisTextBold = 0
slice1Display.PolarAxes.LastRadialAxisTextItalic = 0
slice1Display.PolarAxes.LastRadialAxisTextShadow = 0
slice1Display.PolarAxes.LastRadialAxisTextFontSize = 12
slice1Display.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
slice1Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
slice1Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
slice1Display.PolarAxes.SecondaryRadialAxesTextBold = 0
slice1Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
slice1Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
slice1Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
slice1Display.PolarAxes.EnableDistanceLOD = 1
slice1Display.PolarAxes.DistanceLODThreshold = 0.7
slice1Display.PolarAxes.EnableViewAngleLOD = 1
slice1Display.PolarAxes.ViewAngleLODThreshold = 0.7
slice1Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
slice1Display.PolarAxes.PolarTicksVisibility = 1
slice1Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
slice1Display.PolarAxes.TickLocation = 'Both'
slice1Display.PolarAxes.AxisTickVisibility = 1
slice1Display.PolarAxes.AxisMinorTickVisibility = 0
slice1Display.PolarAxes.ArcTickVisibility = 1
slice1Display.PolarAxes.ArcMinorTickVisibility = 0
slice1Display.PolarAxes.DeltaAngleMajor = 10.0
slice1Display.PolarAxes.DeltaAngleMinor = 5.0
slice1Display.PolarAxes.PolarAxisMajorTickSize = 0.0
slice1Display.PolarAxes.PolarAxisTickRatioSize = 0.3
slice1Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
slice1Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
slice1Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
slice1Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
slice1Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
slice1Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
slice1Display.PolarAxes.ArcMajorTickSize = 0.0
slice1Display.PolarAxes.ArcTickRatioSize = 0.3
slice1Display.PolarAxes.ArcMajorTickThickness = 1.0
slice1Display.PolarAxes.ArcTickRatioThickness = 0.5
slice1Display.PolarAxes.Use2DMode = 0
slice1Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(a5foam, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Plot Data'
plotData1 = PlotData(Input=slice1)

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.UseCache = 0
lineChartView1.ViewSize = [439, 548]
lineChartView1.ChartTitle = ''
lineChartView1.ChartTitleAlignment = 'Center'
lineChartView1.ChartTitleFontFamily = 'Arial'
lineChartView1.ChartTitleFontSize = 14
lineChartView1.ChartTitleBold = 0
lineChartView1.ChartTitleItalic = 0
lineChartView1.ChartTitleColor = [0.0, 0.0, 0.0]
lineChartView1.ShowLegend = 1
lineChartView1.LegendLocation = 'TopRight'
lineChartView1.SortByXAxis = 0
lineChartView1.LegendPosition = [0, 0]
lineChartView1.LegendFontFamily = 'Arial'
lineChartView1.LegendFontSize = 12
lineChartView1.LegendBold = 0
lineChartView1.LegendItalic = 0
lineChartView1.TooltipNotation = 'Mixed'
lineChartView1.TooltipPrecision = 6
lineChartView1.UseOffscreenRendering = 0
lineChartView1.HideTimeMarker = 0
lineChartView1.LeftAxisTitle = ''
lineChartView1.ShowLeftAxisGrid = 1
lineChartView1.LeftAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.LeftAxisColor = [0.0, 0.0, 0.0]
lineChartView1.LeftAxisTitleFontFamily = 'Arial'
lineChartView1.LeftAxisTitleFontSize = 12
lineChartView1.LeftAxisTitleBold = 1
lineChartView1.LeftAxisTitleItalic = 0
lineChartView1.LeftAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.LeftAxisLogScale = 0
lineChartView1.LeftAxisUseCustomRange = 0
lineChartView1.LeftAxisRangeMinimum = 0.0
lineChartView1.LeftAxisRangeMaximum = 1.0
lineChartView1.ShowLeftAxisLabels = 1
lineChartView1.LeftAxisLabelNotation = 'Mixed'
lineChartView1.LeftAxisLabelPrecision = 2
lineChartView1.LeftAxisUseCustomLabels = 0
lineChartView1.LeftAxisLabels = []
lineChartView1.LeftAxisLabelFontFamily = 'Arial'
lineChartView1.LeftAxisLabelFontSize = 12
lineChartView1.LeftAxisLabelBold = 0
lineChartView1.LeftAxisLabelItalic = 0
lineChartView1.LeftAxisLabelColor = [0.0, 0.0, 0.0]
lineChartView1.BottomAxisTitle = ''
lineChartView1.ShowBottomAxisGrid = 1
lineChartView1.BottomAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.BottomAxisColor = [0.0, 0.0, 0.0]
lineChartView1.BottomAxisTitleFontFamily = 'Arial'
lineChartView1.BottomAxisTitleFontSize = 12
lineChartView1.BottomAxisTitleBold = 1
lineChartView1.BottomAxisTitleItalic = 0
lineChartView1.BottomAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.BottomAxisLogScale = 0
lineChartView1.BottomAxisUseCustomRange = 0
lineChartView1.BottomAxisRangeMinimum = 0.0
lineChartView1.BottomAxisRangeMaximum = 1.0
lineChartView1.ShowBottomAxisLabels = 1
lineChartView1.BottomAxisLabelNotation = 'Mixed'
lineChartView1.BottomAxisLabelPrecision = 2
lineChartView1.BottomAxisUseCustomLabels = 0
lineChartView1.BottomAxisLabels = []
lineChartView1.BottomAxisLabelFontFamily = 'Arial'
lineChartView1.BottomAxisLabelFontSize = 12
lineChartView1.BottomAxisLabelBold = 0
lineChartView1.BottomAxisLabelItalic = 0
lineChartView1.BottomAxisLabelColor = [0.0, 0.0, 0.0]
lineChartView1.RightAxisTitle = ''
lineChartView1.ShowRightAxisGrid = 1
lineChartView1.RightAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.RightAxisColor = [0.0, 0.0, 0.0]
lineChartView1.RightAxisTitleFontFamily = 'Arial'
lineChartView1.RightAxisTitleFontSize = 12
lineChartView1.RightAxisTitleBold = 1
lineChartView1.RightAxisTitleItalic = 0
lineChartView1.RightAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.RightAxisLogScale = 0
lineChartView1.RightAxisUseCustomRange = 0
lineChartView1.RightAxisRangeMinimum = 0.0
lineChartView1.RightAxisRangeMaximum = 1.0
lineChartView1.ShowRightAxisLabels = 1
lineChartView1.RightAxisLabelNotation = 'Mixed'
lineChartView1.RightAxisLabelPrecision = 2
lineChartView1.RightAxisUseCustomLabels = 0
lineChartView1.RightAxisLabels = []
lineChartView1.RightAxisLabelFontFamily = 'Arial'
lineChartView1.RightAxisLabelFontSize = 12
lineChartView1.RightAxisLabelBold = 0
lineChartView1.RightAxisLabelItalic = 0
lineChartView1.RightAxisLabelColor = [0.0, 0.0, 0.0]
lineChartView1.TopAxisTitle = ''
lineChartView1.ShowTopAxisGrid = 1
lineChartView1.TopAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.TopAxisColor = [0.0, 0.0, 0.0]
lineChartView1.TopAxisTitleFontFamily = 'Arial'
lineChartView1.TopAxisTitleFontSize = 12
lineChartView1.TopAxisTitleBold = 1
lineChartView1.TopAxisTitleItalic = 0
lineChartView1.TopAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.TopAxisLogScale = 0
lineChartView1.TopAxisUseCustomRange = 0
lineChartView1.TopAxisRangeMinimum = 0.0
lineChartView1.TopAxisRangeMaximum = 1.0
lineChartView1.ShowTopAxisLabels = 1
lineChartView1.TopAxisLabelNotation = 'Mixed'
lineChartView1.TopAxisLabelPrecision = 2
lineChartView1.TopAxisUseCustomLabels = 0
lineChartView1.TopAxisLabels = []
lineChartView1.TopAxisLabelFontFamily = 'Arial'
lineChartView1.TopAxisLabelFontSize = 12
lineChartView1.TopAxisLabelBold = 0
lineChartView1.TopAxisLabelItalic = 0
lineChartView1.TopAxisLabelColor = [0.0, 0.0, 0.0]

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(2, lineChartView1)

# show data in view
plotData1Display = Show(plotData1, lineChartView1)
# trace defaults for the display properties.
plotData1Display.CompositeDataSetIndex = [2]
plotData1Display.AttributeType = 'Point Data'
plotData1Display.UseIndexForXAxis = 1
plotData1Display.XArrayName = 'k'
plotData1Display.SeriesVisibility = ['k (Patches)', 'nut (Patches)', 'omega (Patches)', 'p (Patches)', 'U_Magnitude (Patches)', 'k (airfoil)', 'nut (airfoil)', 'omega (airfoil)', 'p (airfoil)', 'U_Magnitude (airfoil)']
plotData1Display.SeriesLabel = ['k (Patches)', 'k (Patches)', 'nut (Patches)', 'nut (Patches)', 'omega (Patches)', 'omega (Patches)', 'p (Patches)', 'p (Patches)', 'U_X (Patches)', 'U_X (Patches)', 'U_Y (Patches)', 'U_Y (Patches)', 'U_Z (Patches)', 'U_Z (Patches)', 'U_Magnitude (Patches)', 'U_Magnitude (Patches)', 'Points_X (Patches)', 'Points_X (Patches)', 'Points_Y (Patches)', 'Points_Y (Patches)', 'Points_Z (Patches)', 'Points_Z (Patches)', 'Points_Magnitude (Patches)', 'Points_Magnitude (Patches)', 'k (airfoil)', 'k (airfoil)', 'nut (airfoil)', 'nut (airfoil)', 'omega (airfoil)', 'omega (airfoil)', 'p (airfoil)', 'p (airfoil)', 'U_X (airfoil)', 'U_X (airfoil)', 'U_Y (airfoil)', 'U_Y (airfoil)', 'U_Z (airfoil)', 'U_Z (airfoil)', 'U_Magnitude (airfoil)', 'U_Magnitude (airfoil)', 'Points_X (airfoil)', 'Points_X (airfoil)', 'Points_Y (airfoil)', 'Points_Y (airfoil)', 'Points_Z (airfoil)', 'Points_Z (airfoil)', 'Points_Magnitude (airfoil)', 'Points_Magnitude (airfoil)']
plotData1Display.SeriesColor = ['k (Patches)', '0', '0', '0', 'nut (Patches)', '0.89', '0.1', '0.11', 'omega (Patches)', '0.22', '0.49', '0.72', 'p (Patches)', '0.3', '0.69', '0.29', 'U_X (Patches)', '0.6', '0.31', '0.64', 'U_Y (Patches)', '1', '0.5', '0', 'U_Z (Patches)', '0.65', '0.34', '0.16', 'U_Magnitude (Patches)', '0', '0', '0', 'Points_X (Patches)', '0.89', '0.1', '0.11', 'Points_Y (Patches)', '0.22', '0.49', '0.72', 'Points_Z (Patches)', '0.3', '0.69', '0.29', 'Points_Magnitude (Patches)', '0.6', '0.31', '0.64', 'k (airfoil)', '1', '0.5', '0', 'nut (airfoil)', '0.65', '0.34', '0.16', 'omega (airfoil)', '0', '0', '0', 'p (airfoil)', '0.89', '0.1', '0.11', 'U_X (airfoil)', '0.22', '0.49', '0.72', 'U_Y (airfoil)', '0.3', '0.69', '0.29', 'U_Z (airfoil)', '0.6', '0.31', '0.64', 'U_Magnitude (airfoil)', '1', '0.5', '0', 'Points_X (airfoil)', '0.65', '0.34', '0.16', 'Points_Y (airfoil)', '0', '0', '0', 'Points_Z (airfoil)', '0.89', '0.1', '0.11', 'Points_Magnitude (airfoil)', '0.22', '0.49', '0.72']
plotData1Display.SeriesPlotCorner = ['k (Patches)', '0', 'nut (Patches)', '0', 'omega (Patches)', '0', 'p (Patches)', '0', 'U_X (Patches)', '0', 'U_Y (Patches)', '0', 'U_Z (Patches)', '0', 'U_Magnitude (Patches)', '0', 'Points_X (Patches)', '0', 'Points_Y (Patches)', '0', 'Points_Z (Patches)', '0', 'Points_Magnitude (Patches)', '0', 'k (airfoil)', '0', 'nut (airfoil)', '0', 'omega (airfoil)', '0', 'p (airfoil)', '0', 'U_X (airfoil)', '0', 'U_Y (airfoil)', '0', 'U_Z (airfoil)', '0', 'U_Magnitude (airfoil)', '0', 'Points_X (airfoil)', '0', 'Points_Y (airfoil)', '0', 'Points_Z (airfoil)', '0', 'Points_Magnitude (airfoil)', '0']
plotData1Display.SeriesLabelPrefix = ''
plotData1Display.SeriesLineStyle = ['k (Patches)', '1', 'nut (Patches)', '1', 'omega (Patches)', '1', 'p (Patches)', '1', 'U_X (Patches)', '1', 'U_Y (Patches)', '1', 'U_Z (Patches)', '1', 'U_Magnitude (Patches)', '1', 'Points_X (Patches)', '1', 'Points_Y (Patches)', '1', 'Points_Z (Patches)', '1', 'Points_Magnitude (Patches)', '1', 'k (airfoil)', '1', 'nut (airfoil)', '1', 'omega (airfoil)', '1', 'p (airfoil)', '1', 'U_X (airfoil)', '1', 'U_Y (airfoil)', '1', 'U_Z (airfoil)', '1', 'U_Magnitude (airfoil)', '1', 'Points_X (airfoil)', '1', 'Points_Y (airfoil)', '1', 'Points_Z (airfoil)', '1', 'Points_Magnitude (airfoil)', '1']
plotData1Display.SeriesLineThickness = ['k (Patches)', '2', 'nut (Patches)', '2', 'omega (Patches)', '2', 'p (Patches)', '2', 'U_X (Patches)', '2', 'U_Y (Patches)', '2', 'U_Z (Patches)', '2', 'U_Magnitude (Patches)', '2', 'Points_X (Patches)', '2', 'Points_Y (Patches)', '2', 'Points_Z (Patches)', '2', 'Points_Magnitude (Patches)', '2', 'k (airfoil)', '2', 'nut (airfoil)', '2', 'omega (airfoil)', '2', 'p (airfoil)', '2', 'U_X (airfoil)', '2', 'U_Y (airfoil)', '2', 'U_Z (airfoil)', '2', 'U_Magnitude (airfoil)', '2', 'Points_X (airfoil)', '2', 'Points_Y (airfoil)', '2', 'Points_Z (airfoil)', '2', 'Points_Magnitude (airfoil)', '2']
plotData1Display.SeriesMarkerStyle = ['k (Patches)', '0', 'nut (Patches)', '0', 'omega (Patches)', '0', 'p (Patches)', '0', 'U_X (Patches)', '0', 'U_Y (Patches)', '0', 'U_Z (Patches)', '0', 'U_Magnitude (Patches)', '0', 'Points_X (Patches)', '0', 'Points_Y (Patches)', '0', 'Points_Z (Patches)', '0', 'Points_Magnitude (Patches)', '0', 'k (airfoil)', '0', 'nut (airfoil)', '0', 'omega (airfoil)', '0', 'p (airfoil)', '0', 'U_X (airfoil)', '0', 'U_Y (airfoil)', '0', 'U_Z (airfoil)', '0', 'U_Magnitude (airfoil)', '0', 'Points_X (airfoil)', '0', 'Points_Y (airfoil)', '0', 'Points_Z (airfoil)', '0', 'Points_Magnitude (airfoil)', '0']

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotData1Display
plotData1Display.UseIndexForXAxis = 0

# Properties modified on plotData1Display
plotData1Display.XArrayName = 'Points_X'

# Properties modified on plotData1Display
plotData1Display.SeriesVisibility = ['k (Patches)', 'nut (Patches)', 'omega (Patches)', 'p (Patches)', 'U_Magnitude (Patches)', 'k (airfoil)', 'nut (airfoil)', 'p (airfoil)', 'U_Magnitude (airfoil)']
plotData1Display.SeriesColor = ['k (Patches)', '0', '0', '0', 'nut (Patches)', '0.889998', '0.100008', '0.110002', 'omega (Patches)', '0.220005', '0.489998', '0.719997', 'p (Patches)', '0.300008', '0.689998', '0.289998', 'U_X (Patches)', '0.6', '0.310002', '0.639994', 'U_Y (Patches)', '1', '0.500008', '0', 'U_Z (Patches)', '0.650004', '0.340002', '0.160006', 'U_Magnitude (Patches)', '0', '0', '0', 'Points_X (Patches)', '0.889998', '0.100008', '0.110002', 'Points_Y (Patches)', '0.220005', '0.489998', '0.719997', 'Points_Z (Patches)', '0.300008', '0.689998', '0.289998', 'Points_Magnitude (Patches)', '0.6', '0.310002', '0.639994', 'k (airfoil)', '1', '0.500008', '0', 'nut (airfoil)', '0.650004', '0.340002', '0.160006', 'omega (airfoil)', '0', '0', '0', 'p (airfoil)', '0.889998', '0.100008', '0.110002', 'U_X (airfoil)', '0.220005', '0.489998', '0.719997', 'U_Y (airfoil)', '0.300008', '0.689998', '0.289998', 'U_Z (airfoil)', '0.6', '0.310002', '0.639994', 'U_Magnitude (airfoil)', '1', '0.500008', '0', 'Points_X (airfoil)', '0.650004', '0.340002', '0.160006', 'Points_Y (airfoil)', '0', '0', '0', 'Points_Z (airfoil)', '0.889998', '0.100008', '0.110002', 'Points_Magnitude (airfoil)', '0.220005', '0.489998', '0.719997']
plotData1Display.SeriesPlotCorner = ['Points_Magnitude (Patches)', '0', 'Points_Magnitude (airfoil)', '0', 'Points_X (Patches)', '0', 'Points_X (airfoil)', '0', 'Points_Y (Patches)', '0', 'Points_Y (airfoil)', '0', 'Points_Z (Patches)', '0', 'Points_Z (airfoil)', '0', 'U_Magnitude (Patches)', '0', 'U_Magnitude (airfoil)', '0', 'U_X (Patches)', '0', 'U_X (airfoil)', '0', 'U_Y (Patches)', '0', 'U_Y (airfoil)', '0', 'U_Z (Patches)', '0', 'U_Z (airfoil)', '0', 'k (Patches)', '0', 'k (airfoil)', '0', 'nut (Patches)', '0', 'nut (airfoil)', '0', 'omega (Patches)', '0', 'omega (airfoil)', '0', 'p (Patches)', '0', 'p (airfoil)', '0']
plotData1Display.SeriesLineStyle = ['Points_Magnitude (Patches)', '1', 'Points_Magnitude (airfoil)', '1', 'Points_X (Patches)', '1', 'Points_X (airfoil)', '1', 'Points_Y (Patches)', '1', 'Points_Y (airfoil)', '1', 'Points_Z (Patches)', '1', 'Points_Z (airfoil)', '1', 'U_Magnitude (Patches)', '1', 'U_Magnitude (airfoil)', '1', 'U_X (Patches)', '1', 'U_X (airfoil)', '1', 'U_Y (Patches)', '1', 'U_Y (airfoil)', '1', 'U_Z (Patches)', '1', 'U_Z (airfoil)', '1', 'k (Patches)', '1', 'k (airfoil)', '1', 'nut (Patches)', '1', 'nut (airfoil)', '1', 'omega (Patches)', '1', 'omega (airfoil)', '1', 'p (Patches)', '1', 'p (airfoil)', '1']
plotData1Display.SeriesLineThickness = ['Points_Magnitude (Patches)', '2', 'Points_Magnitude (airfoil)', '2', 'Points_X (Patches)', '2', 'Points_X (airfoil)', '2', 'Points_Y (Patches)', '2', 'Points_Y (airfoil)', '2', 'Points_Z (Patches)', '2', 'Points_Z (airfoil)', '2', 'U_Magnitude (Patches)', '2', 'U_Magnitude (airfoil)', '2', 'U_X (Patches)', '2', 'U_X (airfoil)', '2', 'U_Y (Patches)', '2', 'U_Y (airfoil)', '2', 'U_Z (Patches)', '2', 'U_Z (airfoil)', '2', 'k (Patches)', '2', 'k (airfoil)', '2', 'nut (Patches)', '2', 'nut (airfoil)', '2', 'omega (Patches)', '2', 'omega (airfoil)', '2', 'p (Patches)', '2', 'p (airfoil)', '2']
plotData1Display.SeriesMarkerStyle = ['Points_Magnitude (Patches)', '0', 'Points_Magnitude (airfoil)', '0', 'Points_X (Patches)', '0', 'Points_X (airfoil)', '0', 'Points_Y (Patches)', '0', 'Points_Y (airfoil)', '0', 'Points_Z (Patches)', '0', 'Points_Z (airfoil)', '0', 'U_Magnitude (Patches)', '0', 'U_Magnitude (airfoil)', '0', 'U_X (Patches)', '0', 'U_X (airfoil)', '0', 'U_Y (Patches)', '0', 'U_Y (airfoil)', '0', 'U_Z (Patches)', '0', 'U_Z (airfoil)', '0', 'k (Patches)', '0', 'k (airfoil)', '0', 'nut (Patches)', '0', 'nut (airfoil)', '0', 'omega (Patches)', '0', 'omega (airfoil)', '0', 'p (Patches)', '0', 'p (airfoil)', '0']

# Properties modified on plotData1Display
plotData1Display.SeriesVisibility = ['k (Patches)', 'nut (Patches)', 'omega (Patches)', 'p (Patches)', 'U_Magnitude (Patches)', 'k (airfoil)', 'p (airfoil)', 'U_Magnitude (airfoil)']

# Properties modified on plotData1Display
plotData1Display.SeriesVisibility = ['k (Patches)', 'nut (Patches)', 'omega (Patches)', 'p (Patches)', 'U_Magnitude (Patches)', 'p (airfoil)', 'U_Magnitude (airfoil)']

# Properties modified on plotData1Display
plotData1Display.SeriesVisibility = ['k (Patches)', 'nut (Patches)', 'omega (Patches)', 'p (Patches)', 'U_Magnitude (Patches)', 'p (airfoil)']

# Properties modified on lineChartView1
lineChartView1.SortByXAxis = 1

# save data
SaveData('TO_BE_SAVED', proxy=plotData1, Precision=5,
    UseScientificNotation=0,
    WriteAllTimeSteps=0,
    FieldAssociation='Points')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.37698453664779663, -1.6061188244789415, -0.03879009932279587]
renderView1.CameraFocalPoint = [0.37698453664779663, -0.009999990463256836, -0.03879009932279587]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.8855292971021455

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
