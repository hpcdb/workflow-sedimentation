
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
      renderView1.ViewSize = [1121, 415]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [8.999999523162842, 1.0, 1.0000000447034836]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [-5.880407611124431, -12.618806955240096, 7.734768445634158]
      renderView1.CameraFocalPoint = [28.67198714354935, 19.004144605062834, -7.903404238100369]
      renderView1.CameraViewUp = [0.21683030117332994, 0.23200627661466575, 0.9482392673289175]
      renderView1.CameraParallelScale = 12.780607864030875
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1121, height=415, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      test_2xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Slice'
      slice1 = Slice(Input=test_2xmf)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Origin = [9.0, 1.0, 1.0]
      slice1.SliceType.Normal = [0.0, 1.0, 0.0]

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
      test_2xmfDisplay.Representation = 'Outline'
      test_2xmfDisplay.ColorArrayName = ['POINTS', 's']
      test_2xmfDisplay.LookupTable = sLUT
      #test_2xmfDisplay.OSPRayScaleArray = 'u'
      #test_2xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      test_2xmfDisplay.SelectOrientationVectors = 'None'
      test_2xmfDisplay.ScaleFactor = 1.8
      test_2xmfDisplay.SelectScaleArray = 'u'
      test_2xmfDisplay.GlyphType = 'Arrow'
      test_2xmfDisplay.ScalarOpacityUnitDistance = 0.4379839638278351
      #test_2xmfDisplay.GaussianRadius = 0.9
      #test_2xmfDisplay.SetScaleArray = ['POINTS', 'u']
      #test_2xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      #test_2xmfDisplay.OpacityArray = ['POINTS', 'u']
      #test_2xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      #test_2xmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      #test_2xmfDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      #test_2xmfDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # show color legend
      test_2xmfDisplay.SetScalarBarVisibility(renderView1, True)

      # show data from slice1
      slice1Display = Show(slice1, renderView1)
      # trace defaults for the display properties.
      slice1Display.ColorArrayName = ['POINTS', 's']
      slice1Display.LookupTable = sLUT
      #slice1Display.OSPRayScaleArray = 'u'
      #slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      slice1Display.SelectOrientationVectors = 'None'
      slice1Display.ScaleFactor = 0.2
      slice1Display.SelectScaleArray = 'u'
      slice1Display.GlyphType = 'Arrow'
      #slice1Display.GaussianRadius = 0.1
      #slice1Display.SetScaleArray = ['POINTS', 'u']
      #slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
      #slice1Display.OpacityArray = ['POINTS', 'u']
      #slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      #slice1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      #slice1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      #slice1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.7329266667366028, 1.0, 0.5, 0.0]

      # show color legend
      slice1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Position = [0.85, 0.52]
      sLUTColorBar.Position2 = [0.12, 0.42999999999999994]
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(slice1)
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
