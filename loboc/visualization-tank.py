
from paraview.simple import *
from paraview import coprocessing

import datetime as dt

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.0.1 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      start=dt.datetime.now()

      # timeStep = datadescription.GetTimeStep()
      # time = datadescription.GetTime()
      # print "[CATALYST] Visualization  - Time step: " + str(timeStep) + " ; Time: " + str(time)

      # # state file generated using paraview version 5.0.1

      # # ----------------------------------------------------------------
      # # setup views used in the visualization
      # # ----------------------------------------------------------------

      # #### disable automatic camera reset on 'Show'
      # paraview.simple._DisableFirstRenderCameraReset()

      # # Create a new 'Render View'
      # renderView1 = CreateView('RenderView')
      # renderView1.ViewSize = [1174, 808]
      # renderView1.AxesGrid = 'GridAxes3DActor'
      # renderView1.CenterOfRotation = [8.602811544537222, 0.21762272163543006, 1.0218672022980198]
      # renderView1.StereoType = 0
      # renderView1.CameraPosition = [8.602811544537222, 0.21762272163543006, -37.41480343537758]
      # renderView1.CameraFocalPoint = [8.602811544537222, 0.21762272163543006, 16.913457782105223]
      # renderView1.CameraParallelScale = 14.061188690389201
      # renderView1.Background = [0.32, 0.34, 0.43]

      # # register the view with coprocessor
      # # and provide it with information such as the filename to use,
      # # how frequently to write the images, etc.
      # coprocessor.RegisterView(renderView1,
      #     filename='image_%t.png', freq=1, fittoscreen=1, magnification=1, width=1174, height=808, cinema={})
      # renderView1.ViewTime = datadescription.GetTime()

      # # ----------------------------------------------------------------
      # # setup the data processing pipelines
      # # ----------------------------------------------------------------

      # # create a new 'XDMF Reader'
      # # create a producer from a simulation input
      # tank_360_00099xmf = coprocessor.CreateProducer(datadescription, 'input')

      # # ----------------------------------------------------------------
      # # setup color maps and opacity mapes used in the visualization
      # # note: the Get..() functions create a new object, if needed
      # # ----------------------------------------------------------------

      # # get color transfer function/color map for 'd'
      # dLUT = GetColorTransferFunction('d')
      # dLUT.RGBPoints = [-5.14028130308186e-07, 1.0, 1.0, 0.988235, 0.0011214554959083998, 1.0, 1.0, 0.988235, 0.027198748325925996, 0.956862745098039, 0.776470588235294, 0.270588235294118, 0.06459773014695881, 0.952941176470588, 0.772549019607843, 0.250980392156863, 0.10199671407893204, 0.94902, 0.741176, 0.219608, 0.1206962102667999, 0.941176, 0.678431, 0.14902, 0.22439340104722597, 0.921569, 0.545098, 0.054902, 0.2345931108883139, 0.929412, 0.607843, 0.094118, 0.2804918669815468, 0.909804, 0.486275, 0.035294, 0.32469066884517794, 0.890196, 0.411765, 0.019608, 0.39268881938541766, 0.8, 0.305882, 0.0, 0.4207380574863853, 0.760784, 0.239216, 0.0, 0.4207380574863853, 0.717647058823529, 0.207843137254902, 0.00784313725490196, 0.4487872955873531, 0.678431, 0.180392, 0.011765, 0.47683653368832085, 0.6, 0.121569, 0.023529, 0.5048857717892885, 0.501961, 0.054902, 0.031373, 0.5329350098902562, 0.4, 0.039216, 0.058824, 0.5609842479912239, 0.301961, 0.047059, 0.090196]
      # dLUT.ColorSpace = 'Lab'
      # dLUT.NanColor = [1.0, 0.0, 0.0]
      # dLUT.ScalarRangeInitialized = 1.0

      # # get opacity transfer function/opacity map for 'd'
      # dPWF = GetOpacityTransferFunction('d')
      # dPWF.Points = [-5.14028130308186e-07, 0.0, 0.5, 0.0, 0.011899162454889286, 0.322368413209915, 0.5, 0.0, 0.04419828513349719, 0.598684191703796, 0.5, 0.0, 0.11049648353818821, 0.743421077728271, 0.5, 0.0, 0.13429584309036136, 0.736842095851898, 0.5, 0.0, 0.2124937086054963, 0.888157904148102, 0.5, 0.0, 0.3433901565892839, 0.921052634716034, 0.5, 0.0, 0.5609842479912239, 0.559210538864136, 0.5, 0.0]
      # dPWF.ScalarRangeInitialized = 1

      # # ----------------------------------------------------------------
      # # setup the visualization in view 'renderView1'
      # # ----------------------------------------------------------------

      # # show data from tank_360_00099xmf
      # tank_360_00099xmfDisplay = Show(tank_360_00099xmf, renderView1)
      # # trace defaults for the display properties.
      # tank_360_00099xmfDisplay.ColorArrayName = ['POINTS', 'd']
      # tank_360_00099xmfDisplay.LookupTable = dLUT
      # tank_360_00099xmfDisplay.GlyphType = 'Arrow'
      # tank_360_00099xmfDisplay.ScalarOpacityUnitDistance = 0.1425704268507018
      # tank_360_00099xmfDisplay.SetScaleArray = ['POINTS', 'u']
      # tank_360_00099xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      # tank_360_00099xmfDisplay.OpacityArray = ['POINTS', 'u']
      # tank_360_00099xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # # show color legend
      # tank_360_00099xmfDisplay.SetScalarBarVisibility(renderView1, True)

      # # setup the color legend parameters for each legend in this view

      # # get color legend/bar for dLUT in view renderView1
      # dLUTColorBar = GetScalarBar(dLUT, renderView1)
      # dLUTColorBar.Position = [0.8696078431372549, 0.02149938042131351]
      # dLUTColorBar.Title = 'd'
      # dLUTColorBar.ComponentTitle = ''

      # end=dt.datetime.now()
      # elapsedTime = (end.microsecond-start.microsecond)/1e6
      # if(elapsedTime < 0.00000):
      #   elapsedTime = 0.00

      # text_file = open("prov/visualization/paraview-" + str(timeStep) + ".prov", "a+")
      # text_file.write("Visualization:ParaView:Run\n      elapsed-time: %.5f seconds.\n" % (elapsedTime))
      # text_file.close()

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [50]}
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

    start=dt.datetime.now()

    timeStep = datadescription.GetTimeStep()
    time = datadescription.GetTime()
    print "[CATALYST] Visualization  - Time step: " + str(timeStep) + " ; Time: " + str(time)

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    end=dt.datetime.now()
    elapsedTime = (end.microsecond-start.microsecond)/1e6
    if(elapsedTime < 0.00000):
      elapsedTime = 0.00

    text_file = open("prov/visualization/paraview-" + str(timeStep) + ".prov", "a+")
    text_file.write("Visualization:ParaView:Run\n      elapsed-time: %.5f seconds.\n" % (elapsedTime))
    text_file.close()

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

