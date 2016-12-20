
try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
import datetime as dt


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      timeStep = datadescription.GetTimeStep()
      time = datadescription.GetTime()
      print "[CATALYST] Time step: " + str(timeStep) + " ; Time: " + str(time)
      start=dt.datetime.now()

      # 3D analysis
      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'Xdmf3ReaderS'
      # create a producer from a simulation input
      output_2_00003xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Plot Over Line'
      plotOverLine3 = PlotOverLine(Input=output_2_00003xmf,
          Source='High Resolution Line Source')
      plotOverLine3.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine3.Source.Point1 = [4.5, 0.0, 0.0]
      plotOverLine3.Source.Point2 = [4.5, 2.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine4 = PlotOverLine(Input=output_2_00003xmf,
          Source='High Resolution Line Source')
      plotOverLine4.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine4.Source.Point1 = [13.5, 0.0, 0.0]
      plotOverLine4.Source.Point2 = [13.5, 2.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine1 = PlotOverLine(Input=output_2_00003xmf,
          Source='High Resolution Line Source')
      plotOverLine1.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine1.Source.Point1 = [0.0, 1.0, 0.0]
      plotOverLine1.Source.Point2 = [18.0, 1.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine2 = PlotOverLine(Input=output_2_00003xmf,
          Source='High Resolution Line Source')
      plotOverLine2.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine2.Source.Point1 = [9.0, 0.0, 0.0]
      plotOverLine2.Source.Point2 = [9.0, 2.0, 0.0]

      # save data
      SaveData('init_ext_line_0_' + str(timeStep) + ".csv", proxy=plotOverLine1, Precision=5,
        UseScientificNotation=0,
        WriteAllTimeSteps=0,
        FieldAssociation='Points')
      SaveData('init_ext_line_1_' + str(timeStep) + ".csv", proxy=plotOverLine3, Precision=5,
        UseScientificNotation=0,
        WriteAllTimeSteps=0,
        FieldAssociation='Points')
      SaveData('init_ext_line_2_' + str(timeStep) + ".csv", proxy=plotOverLine2, Precision=5,
        UseScientificNotation=0,
        WriteAllTimeSteps=0,
        FieldAssociation='Points')
      SaveData('init_ext_line_3_' + str(timeStep) + ".csv", proxy=plotOverLine4, Precision=5,
        UseScientificNotation=0,
        WriteAllTimeSteps=0,
        FieldAssociation='Points')

      end=dt.datetime.now()
      elapsedTime = (end.microsecond-start.microsecond)/1e6

      text_file = open("prov/rde/data-extraction.prov", "w")
      text_file.write("RDE:DataExtraction:Process\n      elapsed-time: %.5f seconds." % (elapsedTime))
      text_file.close()

      start=dt.datetime.now()

      print "[CATALYST] Visualization"
      # state file generated using paraview version 5.2.0-RC1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1159, 819]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesVisibility = 0
      renderView1.CenterOfRotation = [0.6249999999996305, 1.0, 0.9625]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [26.983990481223646, -9.973797189051998, 6.707412894295681]
      renderView1.CameraFocalPoint = [0.09612421420300904, 10.634109810709727, -2.852800743951642]
      renderView1.CameraViewUp = [-0.24943596481836566, 0.12031750491390951, 0.9608878173160603]
      renderView1.CameraParallelScale = 9.1104335791443
      renderView1.Background = [0.32, 0.34, 0.43]
      renderView1.UseGradientBackground = 1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1159, height=819, cinema={})
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
      sLUT.RGBPoints = [0.5, 1.0, 1.0, 0.988235, 0.501, 1.0, 1.0, 0.988235, 0.525, 0.984314, 0.988235, 0.843137, 0.55, 0.988235, 0.988235, 0.741176, 0.575, 0.980392, 0.968627, 0.654902, 0.6, 0.980392, 0.945098, 0.576471, 0.625, 0.968627, 0.905882, 0.486275, 0.65, 0.968627, 0.862745, 0.388235, 0.675, 0.960784, 0.803922, 0.286275, 0.7, 0.94902, 0.741176, 0.219608, 0.725, 0.941176, 0.678431, 0.14902, 0.75, 0.929412, 0.607843, 0.094118, 0.775, 0.921569, 0.545098, 0.054902, 0.8, 0.909804, 0.486275, 0.035294, 0.825, 0.890196, 0.411765, 0.019608, 0.85, 0.8, 0.305882, 0.0, 0.875, 0.760784, 0.239216, 0.0, 0.9, 0.678431, 0.180392, 0.011765, 0.925, 0.6, 0.121569, 0.023529, 0.95, 0.501961, 0.054902, 0.031373, 0.975, 0.4, 0.039216, 0.058824, 1.0, 0.301961, 0.047059, 0.090196]
      sLUT.ColorSpace = 'Lab'
      sLUT.NanColor = [0.25, 0.0, 0.0]
      sLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 's'
      sPWF = GetOpacityTransferFunction('s')
      sPWF.Points = [0.5, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      sPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mergeBlocks1
      mergeBlocks1Display = Show(mergeBlocks1, renderView1)
      # trace defaults for the display properties.
      mergeBlocks1Display.Representation = 'Outline'
      mergeBlocks1Display.ColorArrayName = ['POINTS', '']
      # mergeBlocks1Display.OSPRayScaleArray = 'u'
      # mergeBlocks1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      mergeBlocks1Display.SelectOrientationVectors = 'None'
      mergeBlocks1Display.ScaleFactor = 1.8
      mergeBlocks1Display.SelectScaleArray = 'u'
      mergeBlocks1Display.GlyphType = 'Arrow'
      mergeBlocks1Display.ScalarOpacityUnitDistance = 1.126649736319026
      # mergeBlocks1Display.GaussianRadius = 0.9
      # mergeBlocks1Display.SetScaleArray = ['POINTS', 'u']
      # mergeBlocks1Display.ScaleTransferFunction = 'PiecewiseFunction'
      # mergeBlocks1Display.OpacityArray = ['POINTS', 'u']
      # mergeBlocks1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show data from clip1
      clip1Display = Show(clip1, renderView1)
      # trace defaults for the display properties.
      clip1Display.Representation = 'Volume'
      clip1Display.ColorArrayName = ['POINTS', 's']
      clip1Display.LookupTable = sLUT
      # clip1Display.OSPRayScaleArray = 'u'
      # clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      clip1Display.SelectOrientationVectors = 'None'
      clip1Display.ScaleFactor = 0.2
      clip1Display.SelectScaleArray = 'u'
      clip1Display.GlyphType = 'Arrow'
      clip1Display.ScalarOpacityUnitDistance = 0.3825478473285972
      # clip1Display.GaussianRadius = 0.1
      # clip1Display.SetScaleArray = ['POINTS', 'u']
      # clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
      # clip1Display.OpacityArray = ['POINTS', 'u']
      # clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      clip1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for sLUT in view renderView1
      sLUTColorBar = GetScalarBar(sLUT, renderView1)
      sLUTColorBar.Position = [0.12658316776396938, 0.11954767726161368]
      sLUTColorBar.Position2 = [0.43000000000000016, 0.11999999999999982]
      sLUTColorBar.Orientation = 'Horizontal'
      sLUTColorBar.Title = 's'
      sLUTColorBar.ComponentTitle = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(clip1)
      # ----------------------------------------------------------------
      
      end=dt.datetime.now()
      elapsedTime = (end.microsecond-start.microsecond)/1e6

      text_file = open("prov/visualization/paraview.prov", "w")
      text_file.write("Visualization:ParaView:Run\n      elapsed-time: %.5f seconds." % (elapsedTime))
      text_file.close()

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
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
coprocessor.EnableLiveVisualization(False)


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
