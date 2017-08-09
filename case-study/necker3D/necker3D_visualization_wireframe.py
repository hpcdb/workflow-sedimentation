
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
      renderView1.ViewSize = [1160, 838]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [5.0, 0.5, 0.5]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [11.94891406603272, 9.625290189204721, 11.833306925544996]
      renderView1.CameraFocalPoint = [5.0, 0.5, 0.5]
      renderView1.CameraViewUp = [-0.22625749770608602, 0.8217839953915441, -0.5229518234503931]
      renderView1.CameraParallelScale = 5.04975246918104
      renderView1.Background = [0.0, 0.0, 0.0]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='wireframe_%t.png', freq=1, fittoscreen=0, magnification=1, width=1160, height=838, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      output_1xmf = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.RGBPoints = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0]
      sLUT.ColorSpace = 'HSV'
      sLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from output_1xmf
      output_1xmfDisplay = Show(output_1xmf, renderView1)
      # trace defaults for the display properties.
      output_1xmfDisplay.Representation = 'Wireframe'
      output_1xmfDisplay.ColorArrayName = ['POINTS', 's']
      output_1xmfDisplay.LookupTable = sLUT
      output_1xmfDisplay.OSPRayScaleArray = 'u'
      output_1xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      output_1xmfDisplay.SelectOrientationVectors = 'None'
      output_1xmfDisplay.SelectScaleArray = 'u'
      output_1xmfDisplay.GlyphType = 'Arrow'
      # output_1xmfDisplay.GlyphTableIndexArray = 'u'
      # output_1xmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
      # output_1xmfDisplay.PolarAxes = 'PolarAxesRepresentation'
      # output_1xmfDisplay.ScalarOpacityFunction = sPWF
      # output_1xmfDisplay.ScalarOpacityUnitDistance = 0.585971866836482
      # output_1xmfDisplay.GaussianRadius = 0.5
      # output_1xmfDisplay.SetScaleArray = ['POINTS', 'u']
      # output_1xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      # output_1xmfDisplay.OpacityArray = ['POINTS', 'u']
      # output_1xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      output_1xmfDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      # sLUTColorBar = GetScalarBar(sLUT, renderView1)
      # sLUTColorBar.Orientation = 'Horizontal'
      # sLUTColorBar.WindowLocation = 'AnyLocation'
      # sLUTColorBar.Position = [0.35611366762028, 0.123383396135885]
      # sLUTColorBar.Title = 's'
      # sLUTColorBar.ComponentTitle = ''
      # sLUTColorBar.ScalarBarLength = 0.330000000000002

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(output_1xmf)
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

    timeStep = datadescription.GetTimeStep()
    time = datadescription.GetTime()
    print "[CATALYST] Visualization (Wireframe)  - Time step: " + str(timeStep) + " ; Time: " + str(time)

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    # coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
    #     image_quality=0, padding_amount=imageFileNamePadding)
    coprocessor.WriteImages(datadescription)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
