
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
      timeStep = datadescription.GetTimeStep()
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1, 1, 1]}
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
    print "[CATALYST] Extraction - Time step: " + str(timeStep) + " ; Time: " + str(time)

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

    # create a new 'XDMF Reader'
    # create a producer from a simulation input
    lock2d_1xmf = coprocessor.CreateProducer(datadescription, 'input')

    # create a new 'Plot Over Line'
    plotOverLine2 = PlotOverLine(Input=lock2d_1xmf,
        Source='High Resolution Line Source')

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine2.Source.Point1 = [4.75, 0.0, 0.0]
    plotOverLine2.Source.Point2 = [4.75, 2.0, 0.0]

    # create a new 'Plot Over Line'
    plotOverLine1 = PlotOverLine(Input=lock2d_1xmf,
        Source='High Resolution Line Source')

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine1.Source.Point2 = [20.0, 0.0, 0.0]

    # create a new 'Plot Over Line'
    plotOverLine3 = PlotOverLine(Input=lock2d_1xmf,
        Source='High Resolution Line Source')

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine3.Source.Point1 = [9.25, 0.0, 0.0]
    plotOverLine3.Source.Point2 = [9.25, 2.0, 0.0]

    # create a new 'Plot Over Line'
    plotOverLine4 = PlotOverLine(Input=lock2d_1xmf,
        Source='High Resolution Line Source')

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine4.Source.Point1 = [13.75, 0.0, 0.0]
    plotOverLine4.Source.Point2 = [13.75, 2.0, 0.0]

    # save data
    SaveData('init_ext_line_0_' + str(timeStep) + ".csv", proxy=plotOverLine1, Precision=5,
      UseScientificNotation=0,
      WriteAllTimeSteps=0,
      FieldAssociation='Points')
    SaveData('init_ext_line_1_' + str(timeStep) + ".csv", proxy=plotOverLine2, Precision=5,
      UseScientificNotation=0,
      WriteAllTimeSteps=0,
      FieldAssociation='Points')
    SaveData('init_ext_line_2_' + str(timeStep) + ".csv", proxy=plotOverLine3, Precision=5,
      UseScientificNotation=0,
      WriteAllTimeSteps=0,
      FieldAssociation='Points')
    SaveData('init_ext_line_3_' + str(timeStep) + ".csv", proxy=plotOverLine4, Precision=5,
      UseScientificNotation=0,
      WriteAllTimeSteps=0,
      FieldAssociation='Points')
