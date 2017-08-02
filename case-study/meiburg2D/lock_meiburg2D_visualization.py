
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.4.1-RC2 64 bits

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.4.1-RC2

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1204, 381]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [10.0, 1.0, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [9.96826879361301, 0.622012116551581, 10000.0]
      renderView1.CameraFocalPoint = [9.96826879361301, 0.622012116551581, 0.0]
      renderView1.CameraParallelScale = 3.20220008704513
      renderView1.Background = [0.0, 0.0, 0.0]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1204, height=381, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      lock2d_1xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Plot Over Line'
      plotOverLine3 = PlotOverLine(Input=lock2d_1xmf,
          Source='High Resolution Line Source')
      plotOverLine3.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine3.Source.Point1 = [9.25, 0.0, 0.0]
      plotOverLine3.Source.Point2 = [9.25, 2.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine1 = PlotOverLine(Input=lock2d_1xmf,
          Source='High Resolution Line Source')
      plotOverLine1.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine1.Source.Point2 = [20.0, 0.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine2 = PlotOverLine(Input=lock2d_1xmf,
          Source='High Resolution Line Source')
      plotOverLine2.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine2.Source.Point1 = [4.75, 0.0, 0.0]
      plotOverLine2.Source.Point2 = [4.75, 2.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine4 = PlotOverLine(Input=lock2d_1xmf,
          Source='High Resolution Line Source')
      plotOverLine4.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine4.Source.Point1 = [13.75, 0.0, 0.0]
      plotOverLine4.Source.Point2 = [13.75, 2.0, 0.0]

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.RGBPoints = [-0.228584237141987, 0.0, 0.0, 0.0, 0.0272793302231074, 0.0, 0.0, 1.0, 0.283142897588202, 0.0, 1.0, 1.0, 0.523955666872996, 0.0, 1.0, 0.0, 0.779819234238091, 1.0, 1.0, 0.0, 1.03568280160318, 1.0, 0.0, 0.0, 1.27649557088798, 0.878431372549, 0.0, 1.0]
      sLUT.ColorSpace = 'RGB'
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.Points = [-0.228584237141987, 0.0, 0.5, 0.0, 1.27649557088798, 1.0, 0.5, 0.0]
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from lock2d_1xmf
      lock2d_1xmfDisplay = Show(lock2d_1xmf, renderView1)
      # trace defaults for the display properties.
      lock2d_1xmfDisplay.Representation = 'Surface'
      lock2d_1xmfDisplay.ColorArrayName = ['POINTS', 's']
      lock2d_1xmfDisplay.LookupTable = sLUT
      lock2d_1xmfDisplay.OSPRayScaleArray = 'u'
      lock2d_1xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      lock2d_1xmfDisplay.SelectOrientationVectors = 'None'
      lock2d_1xmfDisplay.ScaleFactor = 2.0
      lock2d_1xmfDisplay.SelectScaleArray = 'u'
      lock2d_1xmfDisplay.GlyphType = 'Arrow'
      lock2d_1xmfDisplay.GlyphTableIndexArray = 'u'
      lock2d_1xmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
      lock2d_1xmfDisplay.PolarAxes = 'PolarAxesRepresentation'
      lock2d_1xmfDisplay.ScalarOpacityFunction = sPWF
      lock2d_1xmfDisplay.ScalarOpacityUnitDistance = 0.5830923807774
      lock2d_1xmfDisplay.GaussianRadius = 1.0
      lock2d_1xmfDisplay.SetScaleArray = ['POINTS', 'u']
      lock2d_1xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      lock2d_1xmfDisplay.OpacityArray = ['POINTS', 'u']
      lock2d_1xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      lock2d_1xmfDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.WindowLocation = 'AnyLocation'
      sLUTColorBar.Position = [0.272128880824162, 0.122994652406417]
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''
      sLUTColorBar.ScalarBarLength = 0.330000000000002

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(lock2d_1xmf)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1, 1, 1]}
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

    timeStep = datadescription.GetTimeStep()
    time = datadescription.GetTime()
    print "[CATALYST] Visualization (Surface - s) - Time step: " + str(timeStep) + " ; Time: " + str(time)

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
