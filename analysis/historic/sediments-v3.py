
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.0.1 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.0.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1627, 813]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [8.99999988079071, 0.000999927520751953, 1.00000000745058]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [8.968048319187025, 10000.001, 1.0]
      renderView1.CameraFocalPoint = [8.968048319187025, 0.00100000001490116, 1.0]
      renderView1.CameraViewUp = [1.0, 0.0, 0.0]
      renderView1.CameraParallelScale = 12.9883582504745
      renderView1.Background = [0.319996948195621, 0.34000152590219, 0.429999237048905]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1627, height=813, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Xdmf3ReaderS'
      # create a producer from a simulation input
      out__1_xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Plot Over Line'
      plotOverLine1 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine1.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine1.Source.Point1 = [0.0, 1.0, 0.0]
      plotOverLine1.Source.Point2 = [18.0, 1.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine2 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine2.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine2.Source.Point1 = [9.0, 0.0, 0.0]
      plotOverLine2.Source.Point2 = [9.0, 2.0, 0.0]

      # create a new 'Slice'
      slice1 = Slice(Input=out__1_xmf)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Origin = [9.0, 0.001, 0.0]
      slice1.SliceType.Normal = [0.0, 1.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine4 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine4.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine4.Source.Point1 = [13.5, 0.0, 0.0]
      plotOverLine4.Source.Point2 = [13.5, 2.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine3 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine3.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine3.Source.Point1 = [4.5, 0.0, 0.0]
      plotOverLine3.Source.Point2 = [4.5, 2.0, 0.0]

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.RGBPoints = [-0.0015562516571071773, 0.231373, 0.298039, 0.752941, 0.5002170601241533, 0.865003, 0.865003, 0.865003, 1.001990371905418, 0.705882, 0.0156863, 0.14902]
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.Points = [-0.0015562516571071773, 0.0, 0.5, 0.0, 1.001990371905418, 1.0, 0.5, 0.0]
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from slice1
      slice1Display = Show(slice1, renderView1)
      # trace defaults for the display properties.
      slice1Display.ColorArrayName = ['POINTS', 's']
      slice1Display.LookupTable = sLUT
      slice1Display.GlyphType = 'Arrow'
      slice1Display.SetScaleArray = ['POINTS', 'u']
      slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
      slice1Display.OpacityArray = ['POINTS', 'u']
      slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      slice1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Position = [0.214575315195753, 0.0617931456548348]
      sLUTColorBar.Position2 = [0.43, 0.12]
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1, 1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
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
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
