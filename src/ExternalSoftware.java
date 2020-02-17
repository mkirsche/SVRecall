
public class ExternalSoftware {
	
	/*
	 * Runs Sniffles, currently with the set of parameters we use for the rest of the pipeline
	 */
	static void runSniffles(String inputFile, String inputVcf, String outFile) throws Exception
	{
		System.out.println("Running Sniffles on " + inputFile);
		String command = Rescreen.SNIFFLES_PATH + " --min_support 2 --max_distance " + Rescreen.SNIFFLES_MAX_DIST + " --max_num_splits 10 "
				+ "--min_length 20 --num_reads_report -1 --min_seq_size 1000 -m " + inputFile + " -v " + outFile;
		if(Rescreen.SNIFFLES_GENOTYPE)
		{
			command += " --Ivcf " + inputVcf;
		}
		Process child = Runtime.getRuntime().exec(command);
	    int exit = child.waitFor();
	    if(exit != 0)
	    {
	    	System.out.println("Command failed: " + command);
	    	System.out.println("Exit code: " + exit);
	    	System.exit(1);
	    }
	    else
	    {
	    	System.out.println("Sniffles ran successfully");
	    }
	}
	
	/*
	 * Extracts a specific region from a BAM file and makes a new BAM file
	 */
	static void runSamtoolsView(String inputFile, String outFile, String region) throws Exception
	{
		System.out.println("Extracting reads from " + inputFile + " " + region);
		String command = Rescreen.SAMTOOLS_PATH + " view -b -h -o " + outFile + " " + inputFile + " " + region;
		Process child = Runtime.getRuntime().exec(command);
	    int exit = child.waitFor();
	    if(exit != 0)
	    {
	    	System.out.println("Command failed: " + command);
	    	System.out.println("Exit code: " + exit);
	    	System.exit(1);
	    }
	    else
	    {
	    	System.out.println("Samtools ran successfully");
	    }
	}
	
	/*
	 * Indexes a BAM file with Samtools
	 */
	static void runSamtoolsIndex(String inputFile) throws Exception
	{
		System.out.println("Indexing " + inputFile);
		String command = Rescreen.SAMTOOLS_PATH + " index " + inputFile;
		Process child = Runtime.getRuntime().exec(command);
	    int exit = child.waitFor();
	    if(exit != 0)
	    {
	    	System.out.println("Command failed: " + command);
	    	System.out.println("Exit code: " + exit);
	    	System.exit(1);
	    }
	    else
	    {
	    	System.out.println("Samtools ran successfully");
	    }
	}
}
