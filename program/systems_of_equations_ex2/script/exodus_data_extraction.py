#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'ExodusIIReader'
oute = ExodusIIReader(FileName=['/home/vitor/Documents/dev/workflow-sedimentation/program/systems_of_equations_ex2/out.e'])
oute.PointVariables = []
oute.SideSetArrayStatus = []

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on oute
oute.PointVariables = ['vel_', 'p']
oute.ElementBlocks = ['Unnamed block ID: 0 Type: QUAD9']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1611, 832]

# show data in view
outeDisplay = Show(oute, renderView1)
# trace defaults for the display properties.
outeDisplay.ColorArrayName = [None, '']
outeDisplay.OSPRayScaleArray = 'GlobalNodeId'
outeDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
outeDisplay.SelectOrientationVectors = 'GlobalNodeId'
outeDisplay.ScaleFactor = 0.1
outeDisplay.SelectScaleArray = 'GlobalNodeId'
outeDisplay.GlyphType = 'Arrow'
outeDisplay.ScalarOpacityUnitDistance = 0.19193831036664846
outeDisplay.GaussianRadius = 0.05
outeDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
outeDisplay.ScaleTransferFunction = 'PiecewiseFunction'
outeDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
outeDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# set scalar coloring
ColorBy(outeDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
outeDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# set scalar coloring
ColorBy(outeDisplay, ('POINTS', 'p'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
outeDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
outeDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 0.0, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1

# create a new 'Slice'
slice1 = Slice(Input=oute)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.0]

# set active source
SetActiveSource(slice1)

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.ColorArrayName = ['POINTS', 'p']
slice1Display.LookupTable = pLUT
slice1Display.OSPRayScaleArray = 'p'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GaussianRadius = 0.05
slice1Display.SetScaleArray = ['POINTS', 'p']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'p']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(oute, renderView1)

# set active source
SetActiveSource(oute)

# show data in view
outeDisplay = Show(oute, renderView1)

# show color bar/color legend
outeDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
slice1Display = Show(slice1, renderView1)

# hide data in view
Hide(oute, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(slice1)

# set active source
SetActiveSource(oute)

# show data in view
outeDisplay = Show(oute, renderView1)

# show color bar/color legend
outeDisplay.SetScalarBarVisibility(renderView1, True)

# set scalar coloring
ColorBy(outeDisplay, ('POINTS', 'vel_'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
outeDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
outeDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vel_'
vel_LUT = GetColorTransferFunction('vel_')
vel_LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 0.0, 0.705882, 0.0156863, 0.14902]
vel_LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'vel_'
vel_PWF = GetOpacityTransferFunction('vel_')
vel_PWF.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]
vel_PWF.ScalarRangeInitialized = 1

# set active source
SetActiveSource(slice1)

# set active source
SetActiveSource(oute)

# set active source
SetActiveSource(slice1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'vel_'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(oute)

# turn off scalar coloring
ColorBy(outeDisplay, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vel_LUT, renderView1)

# set active source
SetActiveSource(slice1)