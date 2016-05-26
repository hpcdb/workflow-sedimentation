package example;

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
            ProvenanceAnalysis.analyzeElapsedTime(logDirectory);
        }
    }
    
}
