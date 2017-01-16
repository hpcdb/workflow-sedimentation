
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
      renderView1.ViewSize = [843, 548]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [9.0, 1.0, 1.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [-6.952177934155311, -13.115707925050087, 7.750935516468179]
      renderView1.CameraFocalPoint = [15.657178373648106, 10.328436290391561, -5.599344535473394]
      renderView1.CameraViewUp = [0.30271277448995376, 0.23416332177220445, 0.9238682345969048]
      renderView1.CameraParallelScale = 9.1104335791443
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=1, magnification=1, width=843, height=548, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XDMF Reader'
      # create a producer from a simulation input
      necker_24xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Merge Blocks'
      mergeBlocks1 = MergeBlocks(Input=necker_24xmf)

      # create a new 'Slice'
      slice1 = Slice(Input=mergeBlocks1)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Origin = [9.0, 1.0, 1.0]
      slice1.SliceType.Normal = [0.0, 1.0, 0.0]

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'u'
      uLUT = GetColorTransferFunction('u')
      uLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 0.0, 0.705882, 0.0156863, 0.14902]
      uLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'u'
      uPWF = GetOpacityTransferFunction('u')
      uPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]
      uPWF.ScalarRangeInitialized = 1

      # get color transfer function/color map for 's'
      sLUT = GetColorTransferFunction('s')
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
      mergeBlocks1Display.Representation = 'Outline'
      mergeBlocks1Display.ColorArrayName = ['POINTS', 'u']
      mergeBlocks1Display.LookupTable = uLUT
      mergeBlocks1Display.SelectOrientationVectors = 'None'
      mergeBlocks1Display.ScaleFactor = 1.8
      mergeBlocks1Display.SelectScaleArray = 'u'
      mergeBlocks1Display.GlyphType = 'Arrow'
      mergeBlocks1Display.ScalarOpacityUnitDistance = 0.21899198191391755
      

      # show data from slice1
      slice1Display = Show(slice1, renderView1)
      # trace defaults for the display properties.
      slice1Display.ColorArrayName = ['POINTS', 's']
      slice1Display.LookupTable = sLUT
      slice1Display.SelectOrientationVectors = 'None'
      slice1Display.ScaleFactor = 1.8
      slice1Display.SelectScaleArray = 'u'
      slice1Display.GlyphType = 'Arrow'
      
      # show color legend
      slice1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Position = [0.2805106888361045, 0.05445155393053011]
      sLUTColorBar.Position2 = [0.43000000000000005, 0.11999999999999994]
      sLUTColorBar.Orientation = 'Horizontal'
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
  freqs = {'input': [1, 1, 1]}
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