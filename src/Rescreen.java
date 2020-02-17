import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.TreeMap;

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
	static String SAMTOOLS_PATH = "";
	static String SNIFFLES_PATH = "";
	static int MAX_DIST = 500;
	static int PADDING = 10000;
	static boolean SNIFFLES_GENOTYPE = false;
	static int SNIFFLES_MAX_DIST = 1000;
	static void usage()
	{
		System.out.println("java -cp src Rescreen <args>");
		System.out.println();
		System.out.println("Required arguments:");
		System.out.println("  merged_vcf=<string>");
		System.out.println("  bam_file_list=<string>");
		System.out.println("  out_file=<string>");
		System.out.println("Optional arguments:");
		System.out.println("  sniffles_path=<string>");
		System.out.println("  samtools_path=<string>");
		System.out.println("  max_dist=<int>");
	}
	static void parseArgs(String[] args)
	{
		for(String s : args)
		{
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
				if(s.endsWith("sniffles_genotype"))
				{
					SNIFFLES_GENOTYPE = true;
				}
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
			else if(param.equals("sniffles_path"))
			{
				SNIFFLES_PATH = val;
			}
			else if(param.equals("samtools_path"))
			{
				SAMTOOLS_PATH = val;
			}
			else if(param.equals("max_dist"))
			{
				MAX_DIST = Integer.parseInt(val);
			}
			else if(param.equals("padding"))
			{
				PADDING = Integer.parseInt(val);
			}
			else if(param.equals("sniffles_max_dist"))
			{
				SNIFFLES_MAX_DIST = Integer.parseInt(val);
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
	static boolean variantExists(int sample, VcfEntry entry, String vcfFile) throws Exception
	{
		Variant toGenotype = VariantInput.fromVcfEntry(entry, 0);
		Scanner input = new Scanner(new FileInputStream(new File(vcfFile)));
		int count = 0;
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			count++;
			VcfEntry cur = VcfEntry.fromLine(line);
			if(cur.getType().equals(entry.getType()) && cur.getStrand().equals(entry.getStrand()))
			{
				Variant candidateMatch = VariantInput.fromVcfEntry(cur, 1);
				double dist = toGenotype.distance(candidateMatch);
				if(dist < MAX_DIST)
				{
					input.close();
					return true;
				}
			}
		}
		input.close();
		System.out.println("Failed to find the variant in sample " + sample + " with " + count + " calls");
		return false;
	}
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		ArrayList<String> bamFiles = PipelineManager.getFilesFromList(BAM_FILE_LIST);
		
		Scanner input = new Scanner(new FileInputStream(new File(MERGED_VCF)));
		PrintWriter out = new PrintWriter(new File(OUT_FILE));
		ArrayList<String> header = new ArrayList<String>();
		SupportVectorTally svt = new SupportVectorTally();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			if(line.startsWith("#"))
			{
				header.add(line);
				continue;
			}
			VcfEntry entry = VcfEntry.fromLine(line);
			
			// Ignore translocations/breakend variants for now
			if(entry.getType().equals("TRA") || entry.getType().equals("BND"))
			{
				continue;
			}
			
			// Ignore non-specific calls
			if(entry.hasInfoField("IS_SPECIFIC") && entry.getInfo("IS_SPECIFIC").equals("0"))
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
						String gtFile = tmpFile + ".togenotype.vcf";
						String snifflesFile = tmpFile + ".vcf";
						String region = entry.getChromosome() + ":" + Math.max(1, entry.getPos() - PADDING) + "-" + (entry.getPos() + PADDING);
						boolean crashed = false;
						try {
							ExternalSoftware.runSamtoolsView(bamFile, tmpFile, region);
							ExternalSoftware.runSamtoolsIndex(tmpFile);
							if(SNIFFLES_GENOTYPE)
							{
								PrintWriter vcfOut = new PrintWriter(new File(gtFile));
								for(String h : header) vcfOut.println(h);
								vcfOut.println(line);
								vcfOut.close();
							}
							ExternalSoftware.runSniffles(tmpFile, gtFile, snifflesFile);
						} catch (Exception e) {
							System.out.println("Error extracting " + bamFile + " region " + region);
							System.out.println(e.getMessage());
							crashed = true;
						}
						
						// Delete temporary files
						File f;
						if((f = new File(tmpFile)).exists())
						{
							f.delete();
						}
						if((f = new File(gtFile)).exists())
						{
							f.delete();
						}
						if((f = new File(tmpFile + ".bai")).exists())
						{
							f.delete();
						}
						if(crashed)
						{
							continue;
						}
						try {
							if(variantExists(i, entry, snifflesFile))
							{
								System.out.println("Found the variant in sample " + i);
								newSuppVec[i] = '1';
							}
						} catch(Exception e) { 
							System.out.println("Error genotyping " + entry.getId());
							System.out.println(e.getMessage());
						}
						if((f = new File(snifflesFile)).exists())
						{
							f.delete();
						}
						if((f = new File(snifflesFile + "_tmp_genotype")).exists())
						{
							f.delete();
						}
							
					}
				}
			}
			svt.add(new String(newSuppVec));
			if(!(new String(newSuppVec)).equals(suppVec))
			{
				System.out.println(entry.getId()+" "+suppVec+" "+new String(newSuppVec));
				entry.setInfo("SUPP_VEC", new String(newSuppVec));
				out.println(entry);
			}
		}
		System.out.println(svt);
		input.close();
		out.close();
		
	}
	
	static class SupportVectorTally
	{
		TreeMap<String, Integer> counts;
		int size;
		SupportVectorTally()
		{
			size = 0;
			counts = new TreeMap<String, Integer>();
		}
		void add(String suppVec)
		{
			size++;
			counts.putIfAbsent(suppVec, 0);
			counts.put(suppVec, 1 + counts.get(suppVec));
		}
		public String toString()
		{
			StringBuilder res = new StringBuilder("");
			res.append("Number of variants: " + size);
			res.append("Support Vector Counts:");
			for(String s : counts.keySet())
			{
				res.append("  " + s + ": " + counts.get(s));
			}
			return res.toString();
		}
	}
}
