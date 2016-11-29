
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.2.0-RC1 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.2.0-RC1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1065, 819]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesVisibility = 0
      renderView1.CenterOfRotation = [0.6249999999996305, 1.0, 0.9625]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [25.784851348984557, -11.981600320770003, 7.795110138702037]
      renderView1.CameraFocalPoint = [1.004194080245864, 10.664518224366308, -2.7932423432402795]
      renderView1.CameraViewUp = [-0.282945980569641, 0.13441008288664114, 0.9496712598040917]
      renderView1.CameraParallelScale = 9.1104335791443
      renderView1.Background = [0.32, 0.34, 0.43]
      renderView1.UseGradientBackground = 1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=10, fittoscreen=0, magnification=1, width=1065, height=819, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      output_24xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Merge Blocks'
      mergeBlocks1 = MergeBlocks(Input=output_24xmf)

      # create a new 'Clip'
      clip1 = Clip(Input=mergeBlocks1)
      clip1.ClipType = 'Scalar'
      clip1.Scalars = ['POINTS', 's']
      clip1.Value = 0.5

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
      sLUT.RGBPoints = [0.0, 0.0, 0.0, 1.90933e-05, 0.06323519531249999, 0.128363, 0.0636265, 0.0, 0.12647039062499998, 0.193795, 0.111057, 0.0, 0.18970608984375, 0.25976, 0.15987, 0.0, 0.25294128515625, 0.328546, 0.210589, 0.0, 0.31617648046874997, 0.399726, 0.26332, 0.0, 0.37941167578125, 0.472969, 0.318261, 0.0, 0.44264687109375, 0.546245, 0.375827, 0.0, 0.5058823536328124, 0.61745, 0.436719, 0.0, 0.569117765625, 0.685545, 0.501113, 0.0, 0.6323529609374999, 0.749578, 0.568799, 0.0, 0.69558815625, 0.80962, 0.6394, 0.0, 0.7588233515625001, 0.865572, 0.712699, 0.10257, 0.82205905078125, 0.917709, 0.787569, 0.233665, 0.8852942460937501, 0.966914, 0.863138, 0.369608, 0.9485294414062501, 1.0, 0.939405, 0.496104, 1.0078125, 0.999225, 1.0, 0.612275]
      sLUT.ColorSpace = 'Lab'
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0078125, 1.0, 0.5, 0.0]
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mergeBlocks1
      mergeBlocks1Display = Show(mergeBlocks1, renderView1)
      # trace defaults for the display properties.
      mergeBlocks1Display.Representation = 'Outline'
      mergeBlocks1Display.ColorArrayName = ['POINTS', '']
      mergeBlocks1Display.OSPRayScaleArray = 'u'
      mergeBlocks1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      mergeBlocks1Display.SelectOrientationVectors = 'None'
      mergeBlocks1Display.ScaleFactor = 1.8
      mergeBlocks1Display.SelectScaleArray = 'u'
      mergeBlocks1Display.GlyphType = 'Arrow'
      mergeBlocks1Display.ScalarOpacityUnitDistance = 1.126649736319026
      mergeBlocks1Display.GaussianRadius = 0.9
      mergeBlocks1Display.SetScaleArray = ['POINTS', 'u']
      mergeBlocks1Display.ScaleTransferFunction = 'PiecewiseFunction'
      mergeBlocks1Display.OpacityArray = ['POINTS', 'u']
      mergeBlocks1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show data from clip1
      clip1Display = Show(clip1, renderView1)
      # trace defaults for the display properties.
      clip1Display.Representation = 'Volume'
      clip1Display.ColorArrayName = ['POINTS', 's']
      clip1Display.LookupTable = sLUT
      clip1Display.OSPRayScaleArray = 'u'
      clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      clip1Display.SelectOrientationVectors = 'None'
      clip1Display.ScaleFactor = 0.2
      clip1Display.SelectScaleArray = 'u'
      clip1Display.GlyphType = 'Arrow'
      clip1Display.ScalarOpacityUnitDistance = 0.3825478473285972
      clip1Display.GaussianRadius = 0.1
      clip1Display.SetScaleArray = ['POINTS', 'u']
      clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
      clip1Display.OpacityArray = ['POINTS', 'u']
      clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      clip1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Position = [0.2181203007518796, 0.04008557457212712]
      sLUTColorBar.Position2 = [0.43000000000000027, 0.11999999999999984]
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(clip1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [10, 10, 10]}
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
    coprocessor.WriteImages(datadescription, rescale_lookuptable=True)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
