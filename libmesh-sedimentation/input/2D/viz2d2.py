
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.2.0 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.2.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1515, 810]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [9.0, 1.0, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [9.016825267340762, -0.27872031789815965, 34.98732148768454]
      renderView1.CameraFocalPoint = [9.016825267340762, -0.27872031789815965, 0.0]
      renderView1.CameraParallelScale = 5.159122082921754
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='necker2d_%t.png', freq=1, fittoscreen=1, magnification=1, width=1515, height=810, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      test_2xmf = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.NumberOfTableValues = 59
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from test_2xmf
      test_2xmfDisplay = Show(test_2xmf, renderView1)
      # trace defaults for the display properties.
      test_2xmfDisplay.Representation = 'Surface With Edges'
      test_2xmfDisplay.ColorArrayName = ['POINTS', 's']
      test_2xmfDisplay.LookupTable = sLUT
      #test_2xmfDisplay.OSPRayScaleArray = 's'
      #test_2xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      test_2xmfDisplay.SelectOrientationVectors = 'None'
      test_2xmfDisplay.ScaleFactor = 1.8
      test_2xmfDisplay.SelectScaleArray = 'None'
      test_2xmfDisplay.GlyphType = 'Arrow'
      test_2xmfDisplay.ScalarOpacityUnitDistance = 0.7444166721795731


      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      #test_2xmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      #test_2xmfDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      #test_2xmfDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # show color legend
      test_2xmfDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Position = [0.264642857142857, 0.23954106280193238]
      sLUTColorBar.Position2 = [0.43, 0.12000000000000022]
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(test_2xmf)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
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
