package example;

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
public class ProvenanceAnalysis {

    static void analyzeElapsedTime(String logDir) {
        File folder = new File(logDir);
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
                    Logger.getLogger(ProvenanceAnalysis.class.getName()).log(Level.SEVERE, null, ex);
                } catch (IOException ex) {
                    Logger.getLogger(ProvenanceAnalysis.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
        DecimalFormat df2 = new DecimalFormat("###.##");
        System.out.println("Total provenance elapsed time: " + Double.valueOf(df2.format(totalElapsedTime)));
    }

}
