
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.0.1 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.0.1

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'Xdmf3ReaderS'
      # create a producer from a simulation input
      out__1_xmf = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Plot Over Line'
      plotOverLine3 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine3.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine3.Source.Point1 = [4.5, 0.0, 0.0]
      plotOverLine3.Source.Point2 = [4.5, 2.0, 0.0]

      # create a new 'Slice'
      slice1 = Slice(Input=out__1_xmf)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Origin = [9.0, 0.001, 0.0]
      slice1.SliceType.Normal = [0.0, 1.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine2 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine2.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine2.Source.Point1 = [9.0, 0.0, 0.0]
      plotOverLine2.Source.Point2 = [9.0, 2.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine1 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine1.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine1.Source.Point1 = [0.0, 1.0, 0.0]
      plotOverLine1.Source.Point2 = [18.0, 1.0, 0.0]

      # create a new 'Plot Over Line'
      plotOverLine4 = PlotOverLine(Input=out__1_xmf,
          Source='High Resolution Line Source')
      plotOverLine4.Tolerance = 2.22044604925031e-16

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine4.Source.Point1 = [13.5, 0.0, 0.0]
      plotOverLine4.Source.Point2 = [13.5, 2.0, 0.0]
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': []}
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
