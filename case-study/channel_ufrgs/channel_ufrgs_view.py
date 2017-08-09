
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.4.0 64 bits

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=4
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.4.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1588, 474]
      renderView1.AnnotationColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
      renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
      renderView1.CenterOfRotation = [72.5, 0.0, -2.5]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [72.5, -83.9551454345303, -2.5]
      renderView1.CameraFocalPoint = [72.5, 197.76420511448637, -2.5]
      renderView1.CameraViewUp = [0.0, 0.0, 1.0]
      renderView1.CameraParallelScale = 72.9143332959988
      renderView1.Background = [1.0, 1.0, 1.0]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1588, height=474, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      canal_240_00050xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Ghost Cells Generator'
      ghostCellsGenerator1 = GhostCellsGenerator(Input=canal_240_00050xmf)

      # create a new 'Merge Blocks'
      mergeBlocks1 = MergeBlocks(Input=ghostCellsGenerator1)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.LockDataRange = 1
      sLUT.RGBPoints = [0.0, 1.0, 1.0, 0.988235, 0.0019999999999999974, 1.0, 1.0, 0.988235, 0.05, 0.984314, 0.988235, 0.843137, 0.1, 0.988235, 0.988235, 0.741176, 0.15, 0.980392, 0.968627, 0.654902, 0.2, 0.980392, 0.945098, 0.576471, 0.25, 0.968627, 0.905882, 0.486275, 0.3, 0.968627, 0.862745, 0.388235, 0.35, 0.960784, 0.803922, 0.286275, 0.4, 0.94902, 0.741176, 0.219608, 0.45, 0.941176, 0.678431, 0.14902, 0.5, 0.929412, 0.607843, 0.094118, 0.55, 0.921569, 0.545098, 0.054902, 0.6, 0.909804, 0.486275, 0.035294, 0.65, 0.890196, 0.411765, 0.019608, 0.7, 0.8, 0.305882, 0.0, 0.7500000000000001, 0.760784, 0.239216, 0.0, 0.8, 0.678431, 0.180392, 0.011765, 0.8500000000000001, 0.6, 0.121569, 0.023529, 0.9, 0.501961, 0.054902, 0.031373, 0.9500000000000001, 0.4, 0.039216, 0.058824, 1.0, 0.301961, 0.047059, 0.090196]
      sLUT.ColorSpace = 'Lab'
      sLUT.NanColor = [0.25, 0.0, 0.0]
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mergeBlocks1
      mergeBlocks1Display = Show(mergeBlocks1, renderView1)
      # trace defaults for the display properties.
      mergeBlocks1Display.Representation = 'Surface'
      mergeBlocks1Display.AmbientColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.ColorArrayName = ['POINTS', 's']
      mergeBlocks1Display.LookupTable = sLUT
      mergeBlocks1Display.OSPRayScaleArray = 'd'
      mergeBlocks1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      mergeBlocks1Display.SelectOrientationVectors = 'None'
      mergeBlocks1Display.ScaleFactor = 14.5
      mergeBlocks1Display.SelectScaleArray = 'd'
      mergeBlocks1Display.GlyphType = 'Arrow'
      mergeBlocks1Display.GlyphTableIndexArray = 'd'
      mergeBlocks1Display.DataAxesGrid = 'GridAxesRepresentation'
      mergeBlocks1Display.PolarAxes = 'PolarAxesRepresentation'
      mergeBlocks1Display.ScalarOpacityFunction = sPWF
      mergeBlocks1Display.ScalarOpacityUnitDistance = 1.410684695570566
      mergeBlocks1Display.GaussianRadius = 7.25
      mergeBlocks1Display.SetScaleArray = ['POINTS', 'd']
      mergeBlocks1Display.ScaleTransferFunction = 'PiecewiseFunction'
      mergeBlocks1Display.OpacityArray = ['POINTS', 'd']
      mergeBlocks1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      mergeBlocks1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      mergeBlocks1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
      mergeBlocks1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

      # show color legend
      mergeBlocks1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.WindowLocation = 'AnyLocation'
      sLUTColorBar.Position = [0.3463224181360202, 0.10970464135021096]
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''
      sLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
      sLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
      sLUTColorBar.ScalarBarLength = 0.3299999999999998

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(mergeBlocks1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
