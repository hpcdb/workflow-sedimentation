
try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  timeStep = -1;
  time = -1;

  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      timeStep = datadescription.GetTimeStep()
      time = datadescription.GetTime()
      print "[CATALYST] Time step: " + str(timeStep) + " ; Time: " + str(time)

      # 3D analysis
      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'Xdmf3ReaderS'
      # create a producer from a simulation input
      output_2_00000xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Plot Over Line'
      # plotOverLine1 = PlotOverLine(Input=output_2_00000xmf,
      #     Source='High Resolution Line Source')
      # plotOverLine1.Tolerance = 2.22044604925031e-16

      # # init the 'High Resolution Line Source' selected for 'Source'
      # plotOverLine1.Source.Point1 = [0.0, 1.0, 0.0]
      # plotOverLine1.Source.Point2 = [18.0, 1.0, 0.0]

      SaveData('ext_plane_' + str(timeStep) + ".csv", proxy=output_2_00000xmf, Precision=5,
        UseScientificNotation=0,
        WriteAllTimeSteps=0,
        FieldAssociation='Points')

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
