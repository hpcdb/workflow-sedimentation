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
        return totalElapsedTime;
    }
    
    static double indexing(String logDir) {
        File folder = new File(logDir + "/indexing");
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
        return totalElapsedTime;
    }

    static double solver(String logDirectory, double provenanceTime, double rdeTime, double indexingTime) {
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
        
        DecimalFormat df2 = new DecimalFormat("###.#####");
        double solverTime = totalElapsedTime - provenanceTime - rdeTime - indexingTime;
        System.out.println("Solver elapsed time (seconds): " + Double.valueOf(df2.format(solverTime)));
        System.out.println("Provenance elapsed time (seconds): " + Double.valueOf(df2.format(provenanceTime)));
        System.out.println("Raw data extraction elapsed time (seconds): " + Double.valueOf(df2.format(rdeTime)));
        System.out.println("Indexing elapsed time (seconds): " + Double.valueOf(df2.format(indexingTime)));
        System.out.println("Total elapsed time (seconds): " + Double.valueOf(df2.format(totalElapsedTime)));
        
        System.out.println("------------------------------------------------------");
        System.out.println("----------------- Percentual Analysis ----------------");
        System.out.println("------------------------------------------------------");
        double solverPercentual = 100.00 * (solverTime / totalElapsedTime);
        double provenancePercentual = 100.00 * (provenanceTime / totalElapsedTime);
        double rdePercentual = 100.00 * (rdeTime / totalElapsedTime);
        double indexingPercentual = 100.00 * (indexingTime / totalElapsedTime);
        System.out.println("Solver elapsed time (%): " + Double.valueOf(df2.format(solverPercentual)));
        System.out.println("Provenance elapsed time (%): " + Double.valueOf(df2.format(provenancePercentual)));
        System.out.println("RDE elapsed time (%): " + Double.valueOf(df2.format(rdePercentual)));
        System.out.println("Indexing elapsed time (%): " + Double.valueOf(df2.format(indexingPercentual)));
        
        return totalElapsedTime;
    }

}
