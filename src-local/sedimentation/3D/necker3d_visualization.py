
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
      renderView1.ViewSize = [726, 333]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [5.0, 0.5, 0.5]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [9.828922564187891, 3.7279780135085843, 9.38339309809687]
      renderView1.CameraFocalPoint = [5.687865424616647, 0.9598157119731661, 1.765412498178403]
      renderView1.CameraViewUp = [-0.1710739362477735, 0.9522169250711229, -0.25301509034993996]
      renderView1.CameraParallelScale = 5.049752469181039
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=726, height=333, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      output_2xmf = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from output_2xmf
      output_2xmfDisplay = Show(output_2xmf, renderView1)
      # trace defaults for the display properties.
      output_2xmfDisplay.Representation = 'Surface'
      output_2xmfDisplay.ColorArrayName = ['POINTS', 's']
      output_2xmfDisplay.LookupTable = sLUT
      output_2xmfDisplay.OSPRayScaleArray = 'u'
      output_2xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      output_2xmfDisplay.SelectOrientationVectors = 'None'
      output_2xmfDisplay.SelectScaleArray = 'u'
      output_2xmfDisplay.GlyphType = 'Arrow'
      output_2xmfDisplay.GlyphTableIndexArray = 'u'
      output_2xmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
      output_2xmfDisplay.PolarAxes = 'PolarAxesRepresentation'
      output_2xmfDisplay.ScalarOpacityFunction = sPWF
      output_2xmfDisplay.ScalarOpacityUnitDistance = 0.5859718668364818
      output_2xmfDisplay.GaussianRadius = 0.5
      output_2xmfDisplay.SetScaleArray = ['POINTS', 'u']
      output_2xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      output_2xmfDisplay.OpacityArray = ['POINTS', 'u']
      output_2xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      output_2xmfDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.WindowLocation = 'AnyLocation'
      sLUTColorBar.Position = [0.29874701798813597, 0.05133412100625215]
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''
      sLUTColorBar.ScalarBarLength = 0.32999999999999985

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(output_2xmf)
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
