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

    static double provenance(String logDir) {
        File folder = new File(logDir + "/log");
        File[] files = folder.listFiles();
        
        String header = "elapsed-time: ";
        double totalElapsedTime = 0.0;
        
        for (File f : files) {
            if (f.isFile() && f.getName().endsWith(".prov")) {
                try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        if(line.contains(header)){
                            double elapsedTime = Double.valueOf(line.trim()
                                    .replaceAll(header, "").replaceAll(" seconds.", ""));
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
        
        DecimalFormat df2 = new DecimalFormat("###.##");
        System.out.println("Provenance elapsed time (seconds): " + Double.valueOf(df2.format(totalElapsedTime)));
        return totalElapsedTime;
    }
    
    static double rawDataExtraction(String logDir) {
        File folder = new File(logDir + "/rde");
        File[] files = folder.listFiles();
        
        String header = "elapsed-time: ";
        double totalElapsedTime = 0.0;
        
        for (File f : files) {
            if (f.isFile() && f.getName().endsWith(".prov")) {
                try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        if(line.contains(header)){
                            double elapsedTime = Double.valueOf(line.trim()
                                    .replaceAll(header, "").replaceAll(" seconds.", ""));
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
        
        DecimalFormat df2 = new DecimalFormat("###.##");
        System.out.println("Raw data extraction elapsed time (seconds): " + Double.valueOf(df2.format(totalElapsedTime)));
        return totalElapsedTime;
    }

    static double solver(String logDirectory, double provenanceTime, double rdeTime) {
        File folder = new File(logDirectory + "/solver");
        File[] files = folder.listFiles();
        
        String header = "elapsed-time: ";
        double totalElapsedTime = 0.0;
        
        for (File f : files) {
            if (f.isFile() && f.getName().endsWith(".prov")) {
                try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        if(line.contains(header)){
                            double elapsedTime = Double.valueOf(line.trim()
                                    .replaceAll(header, "").replaceAll(" seconds.", ""));
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
        
        DecimalFormat df2 = new DecimalFormat("###.##");
        double solverTime = totalElapsedTime - provenanceTime - rdeTime;
        System.out.println("Solver elapsed time (seconds): " + Double.valueOf(df2.format(solverTime)));
        System.out.println("Total elapsed time (seconds): " + Double.valueOf(df2.format(totalElapsedTime)));
        
        System.out.println("------------------------------------------------------");
        System.out.println("----------------- Percentual Analysis ----------------");
        System.out.println("------------------------------------------------------");
        double provenancePercentual = 100.00 * (provenanceTime / totalElapsedTime);
        double rdePercentual = 100.00 * (rdeTime / totalElapsedTime);
        double solverPercentual = 100.00 * (solverTime / totalElapsedTime);
        System.out.println("provenance elapsed time (%): " + Double.valueOf(df2.format(provenancePercentual)));
        System.out.println("rde elapsed time (%): " + Double.valueOf(df2.format(rdePercentual)));
        System.out.println("solver elapsed time (%): " + Double.valueOf(df2.format(solverPercentual)));
        return totalElapsedTime;
    }

}
