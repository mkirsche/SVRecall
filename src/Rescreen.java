import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

/*
 * Code for rescreening a panel of samples for variants within a particular sample.
 * 
 * It assumes that the SUPP_VEC field from merging has one bit per sample (so no multi-step merging handled yet).
 * It also assumes that the first sample is the one we are screening, so only SV calls unique to that sanple are rechecked.
 */
public class Rescreen {
	static String BAM_FILE_LIST = "";
	static String OUT_FILE = "";
	static String MERGED_VCF = "";
	static void usage()
	{
		System.out.println("java -cp src Rescreen <args>");
		System.out.println();
		System.out.println("Required arguments:");
		System.out.println("  merged_vcf=<string>");
		System.out.println("  vcf_file_list=<string>");
		System.out.println("  bam_file_list=<string>");
		System.out.println("  out_file=<string>");
	}
	static void parseArgs(String[] args)
	{
		for(String s : args)
		{
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
				continue;
			}
			String param = s.substring(0, equalsIdx).toLowerCase();
			String val = s.substring(1 + equalsIdx);
			if(param.equals("merged_vcf"))
			{
				MERGED_VCF = val;
			}
			else if(param.equals("bam_file_list"))
			{
				BAM_FILE_LIST = val;
			}
			else if(param.equals("out_file"))
			{
				OUT_FILE = val;
			}
		}
		if(MERGED_VCF.length() == 0 || BAM_FILE_LIST.length() == 0 || OUT_FILE.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	/*
	 * Very basic genotyping of a variant in a VCF file
	 */
	static boolean variantExists(VcfEntry entry, String vcfFile) throws Exception
	{
		Variant toGenotype = VariantInput.fromVcfEntry(entry, 0);
		Scanner input = new Scanner(new FileInputStream(new File(vcfFile)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith(">"))
			{
				continue;
			}
			VcfEntry cur = VcfEntry.fromLine(vcfFile);
			if(cur.getType().equals(entry.getType()) && cur.getStrand().equals(entry.getStrand()))
			{
				Variant candidateMatch = VariantInput.fromVcfEntry(cur, 1);
				double dist = toGenotype.distance(candidateMatch);
				if(dist < 100)
				{
					input.close();
					return true;
				}
			}
		}
		input.close();
		return false;
	}
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		ArrayList<String> bamFiles = PipelineManager.getFilesFromList(BAM_FILE_LIST);
		
		Scanner input = new Scanner(new FileInputStream(new File(MERGED_VCF)));
		PrintWriter out = new PrintWriter(new File(OUT_FILE));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith(">"))
			{
				continue;
			}
			VcfEntry entry = VcfEntry.fromLine(line);
			
			// Ignore translocations/breakend variants for now
			if(entry.getType().equals("TRA") || entry.getType().equals("BND"))
			{
				continue;
			}
			
			String suppVec = entry.getInfo("SUPP_VEC");
			char[] newSuppVec = suppVec.toCharArray();
			if(suppVec.length() > 1 && suppVec.charAt(0) == '1')
			{
				for(int i = 1; i<suppVec.length(); i++)
				{
					if(suppVec.charAt(i) == '0')
					{
						String bamFile = bamFiles.get(i);
						String tmpFile = entry.getId() + "_" + StringUtils.fileBaseName(bamFile);
						String snifflesFile = tmpFile + ".vcf";
						String region = entry.getChromosome() + ":" + Math.max(1, entry.getPos() - 10000) + "-" + (entry.getPos() + 10000);
						try {
							ExternalSoftware.runSamtoolsView(bamFile, tmpFile, region);
							ExternalSoftware.runSamtoolsIndex(tmpFile);
							ExternalSoftware.runSniffles(tmpFile, snifflesFile);
						} catch (Exception e) {
							System.out.println("Error extracting " + bamFile + " region " + region);
							continue;
						}
						
						// Delete temporary files
						File f;
						if((f = new File(tmpFile)).exists())
						{
							f.delete();
						}
						if((f = new File(tmpFile + ".bai")).exists())
						{
							f.delete();
						}
						try {
							if(variantExists(entry, snifflesFile))
							{
								newSuppVec[i] = '1';
							}
						} catch(Exception e) {}
						if((f = new File(snifflesFile)).exists())
						{
							f.delete();
						}
							
					}
				}
			}
			if(!(new String(newSuppVec)).equals(suppVec))
			{
				System.out.println(entry.getId()+" "+suppVec+" "+new String(newSuppVec));
				entry.setInfo("SUPP_VEC", new String(newSuppVec));
				out.println(entry);
			}
		}
		input.close();
		out.close();
		
	}
}
