package main;

/**
 *
 * @author vitor
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if(args.length > 0){
            String logDirectory = args[0];
            System.out.println("######################################################");
            System.out.println("#################### Log Analysis ####################");
            System.out.println("######################################################");
            double provenanceTime = Analysis.provenance(logDirectory);
            double rdeTime = Analysis.rawDataExtraction(logDirectory);
            Analysis.solver(logDirectory, provenanceTime, rdeTime);
            System.out.println("######################################################");
        }
    }
    
}
