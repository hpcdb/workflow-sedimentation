
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
      timeStep = datadescription.GetTimeStep()
      # time = datadescription.GetTime()
      # print "[CATALYST] Extraction - Time step: " + str(timeStep) + " ; Time: " + str(time)
      # start=dt.datetime.now()

      # # state file generated using paraview version 5.0.1

      # # ----------------------------------------------------------------
      # # setup the data processing pipelines
      # # ----------------------------------------------------------------

      # #### disable automatic camera reset on 'Show'
      # paraview.simple._DisableFirstRenderCameraReset()

      # # create a new 'XDMF Reader'
      # # create a producer from a simulation input
      # tank_360_00201xmf = coprocessor.CreateProducer(datadescription, 'input')

      # # create a new 'Plot Over Line'
      # plotOverLine3 = PlotOverLine(Input=tank_360_00201xmf,
      #     Source='High Resolution Line Source')
      # plotOverLine3.Tolerance = 2.22044604925031e-16

      # # init the 'High Resolution Line Source' selected for 'Source'
      # plotOverLine3.Source.Point1 = [10.0, -9.13, 0.39]
      # plotOverLine3.Source.Point2 = [10.0, 9.56, 0.39]

      # # create a new 'Plot Over Line'
      # plotOverLine2 = PlotOverLine(Input=tank_360_00201xmf,
      #     Source='High Resolution Line Source')
      # plotOverLine2.Tolerance = 2.22044604925031e-16

      # # init the 'High Resolution Line Source' selected for 'Source'
      # plotOverLine2.Source.Point1 = [-1.86, 0.0, 0.39]
      # plotOverLine2.Source.Point2 = [19.56, 0.0, 0.39]

      # # create a new 'Plot Over Line'
      # plotOverLine1 = PlotOverLine(Input=tank_360_00201xmf,
      #     Source='High Resolution Line Source')
      # plotOverLine1.Tolerance = 2.22044604925031e-16

      # # init the 'High Resolution Line Source' selected for 'Source'
      # plotOverLine1.Source.Point1 = [-1.86, -9.13, 0.39]
      # plotOverLine1.Source.Point2 = [19.1, 9.56, 0.39]

      # # create a new 'Plot Over Line'
      # plotOverLine4 = PlotOverLine(Input=tank_360_00201xmf,
      #     Source='High Resolution Line Source')
      # plotOverLine4.Tolerance = 2.22044604925031e-16

      # # init the 'High Resolution Line Source' selected for 'Source'
      # plotOverLine4.Source.Point1 = [10.0, -9.13, 1.0]
      # plotOverLine4.Source.Point2 = [10.0, 9.56, 1.0]

      # # save data
      # SaveData('init_ext_line_0_' + str(timeStep) + ".csv", proxy=plotOverLine1, Precision=5,
      #   UseScientificNotation=0,
      #   WriteAllTimeSteps=0,
      #   FieldAssociation='Points')
      # SaveData('init_ext_line_1_' + str(timeStep) + ".csv", proxy=plotOverLine3, Precision=5,
      #   UseScientificNotation=0,
      #   WriteAllTimeSteps=0,
      #   FieldAssociation='Points')
      # SaveData('init_ext_line_2_' + str(timeStep) + ".csv", proxy=plotOverLine2, Precision=5,
      #   UseScientificNotation=0,
      #   WriteAllTimeSteps=0,
      #   FieldAssociation='Points')
      # SaveData('init_ext_line_3_' + str(timeStep) + ".csv", proxy=plotOverLine4, Precision=5,
      #   UseScientificNotation=0,
      #   WriteAllTimeSteps=0,
      #   FieldAssociation='Points')

      # end=dt.datetime.now()
      # elapsedTime = (end.microsecond-start.microsecond)/1e6
      # if(elapsedTime < 0.00000):
      #   elapsedTime = 0.00

      # text_file = open("prov/rde/data-extraction-" + str(timeStep) + ".prov", "a+")
      # text_file.write("RDE:DataExtraction:Process\n      elapsed-time: %.5f seconds.\n" % (elapsedTime))
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

    timeStep = datadescription.GetTimeStep()
    time = datadescription.GetTime()
    print "[CATALYST] Extraction - Time step: " + str(timeStep) + " ; Time: " + str(time)
    start=dt.datetime.now()

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

    # create a new 'XDMF Reader'
    # create a producer from a simulation input
    tank_360_00201xmf = coprocessor.CreateProducer(datadescription, 'input')

    # create a new 'Plot Over Line'
    plotOverLine3 = PlotOverLine(Input=tank_360_00201xmf,
        Source='High Resolution Line Source')
    plotOverLine3.Tolerance = 2.22044604925031e-16

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine3.Source.Point1 = [10.0, -9.13, 0.39]
    plotOverLine3.Source.Point2 = [10.0, 9.56, 0.39]

    # create a new 'Plot Over Line'
    plotOverLine2 = PlotOverLine(Input=tank_360_00201xmf,
        Source='High Resolution Line Source')
    plotOverLine2.Tolerance = 2.22044604925031e-16

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine2.Source.Point1 = [-1.86, 0.0, 0.39]
    plotOverLine2.Source.Point2 = [19.56, 0.0, 0.39]

    # create a new 'Plot Over Line'
    plotOverLine1 = PlotOverLine(Input=tank_360_00201xmf,
        Source='High Resolution Line Source')
    plotOverLine1.Tolerance = 2.22044604925031e-16

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine1.Source.Point1 = [-1.86, -9.13, 0.39]
    plotOverLine1.Source.Point2 = [19.1, 9.56, 0.39]

    # create a new 'Plot Over Line'
    plotOverLine4 = PlotOverLine(Input=tank_360_00201xmf,
        Source='High Resolution Line Source')
    plotOverLine4.Tolerance = 2.22044604925031e-16

    # init the 'High Resolution Line Source' selected for 'Source'
    plotOverLine4.Source.Point1 = [10.0, -9.13, 1.0]
    plotOverLine4.Source.Point2 = [10.0, 9.56, 1.0]

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
    if(elapsedTime < 0.00000):
      elapsedTime = 0.00

    text_file = open("prov/rde/data-extraction-" + str(timeStep) + ".prov", "a+")
    text_file.write("RDE:DataExtraction:Process\n      elapsed-time: %.5f seconds.\n" % (elapsedTime))
    text_file.close()
