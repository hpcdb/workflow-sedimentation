
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.1.0 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.1.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [744, 757]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [8.07329142093658, 1.00000001490116, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [8.07329142093658, 1.00000001490116, 25.8274795258001]
      renderView1.CameraFocalPoint = [8.07329142093658, 1.00000001490116, 0.0]
      renderView1.CameraParallelScale = 11.8422538798836
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=744, height=757, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Xdmf3ReaderS'
      # create a producer from a simulation input
      output_2_00003xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Clip'
      clip1 = Clip(Input=output_2_00003xmf)
      clip1.ClipType = 'Plane'
      clip1.Scalars = ['POINTS', 'u']
      clip1.Value = 1.45877243281148e-06

      # init the 'Plane' selected for 'ClipType'
      clip1.ClipType.Origin = [5.0, 1.0, 0.0]
      clip1.ClipType.Normal = [-1.0, 0.0, 0.0]

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'v'
      vLUT = GetColorTransferFunction('v')
      vLUT.RGBPoints = [-0.0147519137461937, 0.231373, 0.298039, 0.752941, -0.000655078702926521, 0.865003, 0.865003, 0.865003, 0.0134417563403407, 0.705882, 0.0156863, 0.14902]
      vLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'v'
      vPWF = GetOpacityTransferFunction('v')
      vPWF.Points = [-0.0147519137461937, 0.0, 0.5, 0.0, 0.0134417563403407, 1.0, 0.5, 0.0]
      vPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from clip1
      clip1Display = Show(clip1, renderView1)
      # trace defaults for the display properties.
      clip1Display.ColorArrayName = ['POINTS', 'v']
      clip1Display.LookupTable = vLUT
      clip1Display.OSPRayScaleArray = 'u'
      clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      clip1Display.GlyphType = 'Arrow'
      clip1Display.ScalarOpacityUnitDistance = 0.755125263345945
      clip1Display.SetScaleArray = ['POINTS', 'u']
      clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
      clip1Display.OpacityArray = ['POINTS', 'u']
      clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      clip1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for vLUT in view renderView1
      vLUTColorBar = GetScalarBar(vLUT, renderView1)
      vLUTColorBar.Title = 'v'
      vLUTColorBar.ComponentTitle = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(None)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1]}
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
