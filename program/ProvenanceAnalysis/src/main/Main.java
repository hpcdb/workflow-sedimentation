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
            Integer visualizationProcesses = Integer.parseInt(args[1]);
            System.out.println("######################################################");
            System.out.println("#################### Log Analysis ####################");
            System.out.println("######################################################");
            Analysis analysis = new Analysis();
            analysis.setLogDirectory(logDirectory);
            analysis.run(visualizationProcesses);
            analysis.print();
            System.out.println("######################################################");
        }
    }
    
}
