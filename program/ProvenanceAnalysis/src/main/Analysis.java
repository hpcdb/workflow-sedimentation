package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author vitor
 */
public class Analysis {

//    standard
    protected String logDirectory;
    protected double provenanceTime;
    protected double paraviewTime;
    protected double rdeTime;
    protected double rdiTime;
    protected double indexingTime;
    protected double visualizationTime;
//    composed
    double paraviewOverhead;
    double rdiOverhead;
    double solverTime;

    private static double getElapsedTime(String logDir, String subDir) {
        File folder = new File(logDir + "/" + subDir);
        File[] files = folder.listFiles();

        String header = "elapsed-time:";
        double totalElapsedTime = 0.0;

        for (File f : files) {
            if (f.isFile() && f.getName().endsWith(".prov")) {
                try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        if (line.contains(header)) {
                            double elapsedTime = Double.valueOf(line.trim()
                                    .replaceAll(header, "").replaceAll("seconds.", "").trim());
                            totalElapsedTime += elapsedTime;
                        }
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
                } catch (IOException ex) {
                    Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        return totalElapsedTime;
    }
    
    private static double getMaxElapsedTime(String logDir, String subDir) {
        File folder = new File(logDir + "/" + subDir);
        File[] files = folder.listFiles();

        String header = "elapsed-time:";
        double elapsedTime = 0.0;

        for (File f : files) {
            if (f.isFile() && f.getName().endsWith(".prov")) {
                try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                    String line;
                    double maxElapsedTime = 0.0;
                    while ((line = br.readLine()) != null) {
                        if (line.contains(header)) {
                            double currentElapsedTime = Double.valueOf(line.trim()
                                    .replaceAll(header, "").replaceAll("seconds.", "").trim());
                            if(currentElapsedTime > maxElapsedTime){
                                maxElapsedTime = currentElapsedTime;
                            }
                        }
                    }
                    elapsedTime += maxElapsedTime;
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
                } catch (IOException ex) {
                    Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        return elapsedTime;
    }

    private void calculateOverheads() {
        paraviewOverhead = paraviewTime - rdeTime - visualizationTime;
        rdiOverhead = rdiTime - indexingTime;
    }

    private void setProvenanceTime() {
        this.provenanceTime = Analysis.getElapsedTime(logDirectory, "log");
    }

    private void setParaviewTime() {
        this.paraviewTime = Analysis.getMaxElapsedTime(logDirectory, "paraview");
    }

    private void setRDETime() {
        this.rdeTime = Analysis.getMaxElapsedTime(logDirectory, "rde");
    }

    private void setRDITime() {
        this.rdiTime = Analysis.getMaxElapsedTime(logDirectory, "rdi");
    }

    private void setIndexing() {
        this.indexingTime = Analysis.getMaxElapsedTime(logDirectory, "indexing");
    }

    private void setVisualization(Integer visualizationProcesses) {
        this.visualizationTime = visualizationProcesses * Analysis.getMaxElapsedTime(logDirectory, "visualization");
    }

    void setLogDirectory(String logDirectory) {
        this.logDirectory = logDirectory;
    }

    void run(Integer visualizationProcesses) {
        if (logDirectory != null) {
            this.setLogDirectory(logDirectory);
            this.setProvenanceTime();
            this.setParaviewTime();
            this.setRDETime();
            this.setRDITime();
            this.setIndexing();
            this.setVisualization(visualizationProcesses);
            this.calculateOverheads();
            this.calculateSolverTime();
        }
    }

    private void calculateSolverTime() {
        double totalElapsedTime = Analysis.getElapsedTime(logDirectory, "solver");
        solverTime = totalElapsedTime - provenanceTime - paraviewTime - rdiTime;
    }

    void print() {
        if (logDirectory != null) {
            DecimalFormat df2 = new DecimalFormat("###.##");
            double totalElapsedTime = solverTime + provenanceTime
                    + rdeTime + visualizationTime + indexingTime
                    + paraviewOverhead + rdiOverhead;
            System.out.println("Solver elapsed time (seconds): " + Double.valueOf(df2.format(solverTime)));
            System.out.println("Provenance elapsed time (seconds): " + Double.valueOf(df2.format(provenanceTime)));
            System.out.println("Raw data extraction elapsed time (seconds): " + Double.valueOf(df2.format(rdeTime)));
            System.out.println("Visualization elapsed time (seconds): " + Double.valueOf(df2.format(visualizationTime)));
            System.out.println("Raw data indexing elapsed time (seconds): " + Double.valueOf(df2.format(indexingTime)));
            System.out.println("ParaView overhead elapsed time (seconds): " + Double.valueOf(df2.format(paraviewOverhead)));
            System.out.println("RDI overhead  elapsed time (seconds): " + Double.valueOf(df2.format(rdiOverhead)));
            System.out.println("Total elapsed time (seconds): " + Double.valueOf(df2.format(totalElapsedTime)));

            System.out.println("------------------------------------------------------");
            System.out.println("----------------- Percentual Analysis ----------------");
            System.out.println("------------------------------------------------------");
            double solverPercentual = 100.00 * (solverTime / totalElapsedTime);
            double provenancePercentual = 100.00 * (provenanceTime / totalElapsedTime);
            double rdePercentual = 100.00 * (rdeTime / totalElapsedTime);
            double visualizationPercentual = 100.00 * (visualizationTime / totalElapsedTime);
            double indexingPercentual = 100.00 * (indexingTime / totalElapsedTime);
            double paraviewOverheadPercentual = 100.00 * (paraviewOverhead / totalElapsedTime);
            double rdiOverheadPercentual = 100.00 * (rdiOverhead / totalElapsedTime);
            double sumPercentuals = solverPercentual + provenancePercentual + rdePercentual + visualizationPercentual + 
                    indexingPercentual + paraviewOverheadPercentual + rdiOverheadPercentual;
            System.out.println("Solver elapsed time (%): " + Double.valueOf(df2.format(solverPercentual)));
            System.out.println("Provenance elapsed time (%): " + Double.valueOf(df2.format(provenancePercentual)));
            System.out.println("Raw data extraction elapsed time (%): " + Double.valueOf(df2.format(rdePercentual)));
            System.out.println("Visualization elapsed time (%): " + Double.valueOf(df2.format(visualizationPercentual)));
            System.out.println("Raw data indexing elapsed time (%): " + Double.valueOf(df2.format(indexingPercentual)));
            System.out.println("Paraview OVERHEAD elapsed time (%): " + Double.valueOf(df2.format(paraviewOverheadPercentual)));
            System.out.println("Raw data indexing OVERHEAD elapsed time (%): " + Double.valueOf(df2.format(rdiOverheadPercentual)));
            System.out.println("Total elapsed time (%): " + Double.valueOf(df2.format(sumPercentuals)));
        }
    }

}
