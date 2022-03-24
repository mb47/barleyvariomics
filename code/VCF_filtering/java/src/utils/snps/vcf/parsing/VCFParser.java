package utils.snps.vcf.parsing;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import utils.entities.*;

public class VCFParser
{
	public String [] mainHeaders;
	public LinkedList<String> doubleHashHeaders = new LinkedList<String>();
	public static ArrayList<String> sampleNamesList;	
	static boolean debug = false;
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public LinkedList<SNP> parseAndFilterSNPsFromVCF(File vcfFile, int minAlleleOccurrence,float minimumAlternateAlleleFraction, float percentDiffAlleleQuals, float qualScoreCutof)
	{
		//this holds the actual SNP objects
		LinkedList<SNP> snps = new LinkedList<SNP>();
		//these are the VCF entries we need to convert into SNP objects
		System.out.println("parsing VCF file");
		LinkedList<VCFEntry> vcfEntries =  parseVCF(vcfFile);
		
		//counts
		int numMultiAllelicSNPs = 0;
		int numSNPsFailedFilter = 0;
		int count = 0;
		
		for (VCFEntry vcfEntry : vcfEntries)
		{
			count++;
			System.out.print("\rnum SNPs processed: " + count);
			
			SNP snp = configureVCFEntry(vcfEntry, minAlleleOccurrence,minimumAlternateAlleleFraction, numMultiAllelicSNPs);	

			//check the SNP is not null at this point
			//this will happen if it is multi-allelic or a complex haplotype, e.g. CGA
			if(snp == null)
			{
				numMultiAllelicSNPs++;
				continue;
			}
			
			//check the SNP by the criteria passed in as parameters (minAlleleOccurrence, minimumAlternateAlleleFraction)
			boolean pass = filterSNP(snp, minAlleleOccurrence, minimumAlternateAlleleFraction, percentDiffAlleleQuals, qualScoreCutof);	
			
			//if it has passed, add it to the list of SNPs returned
			if(pass)
				snps.add(snp);
			else
				numSNPsFailedFilter++;
		}
	
		System.out.println();
		System.out.println("# SNPs discarded due to failed filters: " + numSNPsFailedFilter);
		System.out.println("# SNPs discarded due to multiple alleles or complex haplotypes: " + numMultiAllelicSNPs);
		
		return snps;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static SNP configureVCFEntry(VCFEntry vcfEntry, int minAlleleOccurrence,float minimumAlternateAlleleFraction, int numMultiAllelicSNPs)
	{
		SNP snp = new SNP();			
		
		//configure this SNP with the info from the VCF entry
		snp.refSeqName = vcfEntry.contig;
		snp.position = vcfEntry.position;
		snp.name = snp.refSeqName + "_" + snp.position;
		snp.qualityScore = vcfEntry.qual;
		
//		System.out.println("\n=========parsing SNP at pos " + snp.position + " in contig " + snp.refSeqName);
		
		//disregard multiallelic SNPS
		if(vcfEntry.ref.length() > 1 || vcfEntry.alt.length() > 1)
		{	
			numMultiAllelicSNPs++;
			return null;
		}
		//assign alleles
		snp.refSeqAllele = new Allele(vcfEntry.ref.charAt(0));
		snp.alternateAllele = new Allele(vcfEntry.alt.charAt(0));
		//determine whether the reference allele or the alternate allele are the minor allele
		snp.alleleList.add(snp.refSeqAllele);
		snp.alleleList.add(snp.alternateAllele);
		
		//parse info string and put into map
		HashMap<String, String> infoMap = new HashMap<String, String>();
		//info string looks like this:
		//AB=0.1;ABP=16.9077;AC=57;AF=0.59375;AN=96;AO=395;CIGAR=1X;DP=561;DPRA=1.68592;EPP=860.742;EPPR=363.475;HWE=-94.3876;LEN=1;MEANALT=1;MQM=255;MQMR=255;NS=48;NUMALT=1;ODDS=0.980829;PAIRED=0;PAIREDR=0;RO=166;RPP=860.742;RPPR=363.475;RUN=1;SAP=860.742;SRP=363.475;TYPE=snp;XAI=0;XAM=0;XAS=0;XRI=0;XRM=0.000847139;XRS=0.000847139;BVAR
		String [] keyValuePairs = vcfEntry.info.split(";");
		if(debug) System.out.println("parsing map");
		for (int i = 0; i < keyValuePairs.length; i++)
		{
			String [] tokens = keyValuePairs[i].split("=");
			//strip trailing ";" off the value
			//e.g. AB=0.1;
			String value;
			//only proceed if we actually had a "=" in the string
			//in some cases we don't but these can be ignored
			if (tokens.length > 1)
			{
				value = null;
//				System.out.println("tokens[1] = " + tokens[1]);
				if (tokens[1].endsWith(";"))
					value = tokens[1].substring(0, (tokens[1].length()));
				else
					value = tokens[1];
				infoMap.put(tokens[0], value);
//				System.out.println("key, value : " + tokens[0] + ", " + value);
			}
		}
		
		//get counts
		snp.refSeqAllele.count = Integer.parseInt(infoMap.get("RO"));
		snp.alternateAllele.count = Integer.parseInt(infoMap.get("AO"));
		//work out which is major and minor allele
		snp.workOutMinorMajorAlleles();
		snp.totalReadCount = Integer.parseInt(infoMap.get("DP"));
		
		//populate list of samples for this SNP
		if(debug) System.out.println("vcfEntry.sampleInfos.size() = " + vcfEntry.sampleInfos.size());		
		processSampleInfos(vcfEntry, snp, minAlleleOccurrence, minimumAlternateAlleleFraction);	
		
		//work out the percentage of filtered heterozygous samples in this SNP
		snp.pctFilteredHets = ((snp.numFilteredHets/(float)snp.samplesList.size())*100);
		
		snp.meanQualAltAllele = snp.altAlleleQualSum/(float)snp.altAlleleCount;
		snp.meanQualRefAllele = snp.refAlleleQualSum/(float)snp.refAlleleCount;
		
		//now compute the statistic we use for filtering on this
		float diffRefQualAltQual = snp.meanQualRefAllele - snp.meanQualAltAllele;
		snp.percentDiffRefAlleleQualAltAlleleQual = (diffRefQualAltQual/snp.meanQualRefAllele)*100;
		
		return snp;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static boolean filterSNP(SNP snp,  int minAlleleOccurrence,float minimumAlternateAlleleFraction, float percentDiffAlleleQualsCutoff, float qualScoreCutoff)
	{
		//check the SNP by the criteria passed in as parameters (minAlleleOccurrence, minimumAlternateAlleleFraction)
		snp.passMAC = snp.minorAllele.count >= minAlleleOccurrence;
		snp.passMAF = snp.getMinorAlleleFrequency() >= minimumAlternateAlleleFraction;
		snp.passPctDiff = snp.percentDiffRefAlleleQualAltAlleleQual <= percentDiffAlleleQualsCutoff;
		snp.passQualScoreFilter = snp.qualityScore >= qualScoreCutoff;
		
		boolean pass = snp.passMAC && snp.passMAF && snp.passPctDiff && snp.passQualScoreFilter;		
		
		if(debug) 
		{
			System.out.println("\n============ filtering =============");	
			System.out.println("snp.minorAllele.count = " + snp.minorAllele.count);	
			System.out.println("minAlleleOccurrence = " + minAlleleOccurrence);	
			System.out.println("snp.getMinorAlleleFrequency() = " + snp.getMinorAlleleFrequency());	
			System.out.println("minimumAlternateAlleleFraction = " + minimumAlternateAlleleFraction);	
			System.out.println("snp.meanQualAltAllele = " + snp.meanQualAltAllele);
			System.out.println("snp.meanQualRefAllele = " + snp.meanQualRefAllele);
			System.out.println("snp.percentDiffRefAlleleQualAltAlleleQual = " + snp.percentDiffRefAlleleQualAltAlleleQual);
			System.out.println("passPctDiff = " + snp.passPctDiff);
			System.out.println("qualScoreCutoff = " + qualScoreCutoff);
			System.out.println("snp.qualityScore = " + snp.qualityScore);
			System.out.println("pass = " + pass);	
			System.out.println("============ ");	
		}

		return pass;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	private static void processSampleInfos(VCFEntry vcfEntry, SNP snp, int minAlleleOccurrence,float minimumAlternateAlleleFraction)
	{
		//populate list of samples for this SNP
		//iterate over all sample infos
		//if the text is not "." then we add the sample to our list
		for(int i =0; i < vcfEntry.sampleInfos.size(); i++)
		{
			String sampleInfo = vcfEntry.sampleInfos.get(i);
			String sampleName =  sampleNamesList.get(i);
			
			if(debug) System.out.println("\n+++++sampleIndex = " + i + ", sample " + sampleName);
			
			if(!sampleInfo.equals("."))
			{
				//get the name of this sample
				SNPSample snpSample = new SNPSample();
				snpSample.name = sampleName;
				//add it to the list
				snp.samplesList.add(snpSample);
				
				//extract sample info
				HashMap<String, String> lookup = getSampleInfoLookup(vcfEntry.format, sampleInfo);
				if(debug) System.out.println("info lookup:");
				for(String key : lookup.keySet())
					System.out.println(key + " : " + lookup.get(key));
				
				//extract alleles and counts and add to list 
				//ref allele
				Allele refAllele = new Allele(snp.refSeqAllele.value);
				refAllele.count = Integer.parseInt(lookup.get("RO"));
				if(refAllele.count > 0)
					snpSample.alleleList.add(refAllele);
				//alt allele
				Allele altAllele = new Allele(snp.alternateAllele.value);
				altAllele.count = Integer.parseInt(lookup.get("AO"));
				if(altAllele.count > 0)
					snpSample.alleleList.add(altAllele);
				
				//work out minor/major of these
				snpSample.workOutMinorMajorAlleles();
				
				//keep track of base qualities at SNP level so we can use this for filtering later
				int refAlleleQualSum = Integer.parseInt(lookup.get("QR"));
				int altAlleleQualSum = Integer.parseInt(lookup.get("QA"));
				
				snp.refAlleleCount += refAllele.count;
				snp.altAlleleCount += altAllele.count;
				snp.refAlleleQualSum += refAlleleQualSum;
				snp.altAlleleQualSum += altAlleleQualSum;
				
				if (debug)
				{
					System.out.println("altAllele.count = " + altAllele.count);
					System.out.println("refAllele.count = " + refAllele.count);
					
					System.out.println("altAlleleQualSum = " + altAlleleQualSum);
					System.out.println("refAlleleQualSum = " + refAlleleQualSum);
					
					System.out.println("majorAllele = " + snpSample.majorAllele + ", count = " + snpSample.majorAllele.count);
					if(snpSample.minorAllele != null) System.out.println("minorAllele = " + snpSample.minorAllele + ", count = " + snpSample.minorAllele.count);
					System.out.println("snpSample.isHomozygous() = " + snpSample.rawGenotypeIsHomozygous());
					System.out.println("sampleInfo = " + sampleInfo);
				}
				
				//ditto for the smallest USABLE allele count in this sample
				//that is, for homozygotes we use the major count, for hets the minor count
				//this also counts the number of hets and homs
				if (snpSample.rawGenotypeIsHomozygous())
					snp.numCleanHomoZygotes++;
				else
				{
					snp.numHets++;

					//this is a value we use for getting a better idea of genuine hets
					//for this we filter the allele counts within a given SNP x sample combination the same way we filter the overall allele frequencies
					//i.e. if we filter SNP below a certain minor allele frequency and minro allele count then we do the same within samples here
					if (snpSample.getMinorAlleleCount() >= minAlleleOccurrence && snpSample.getMinorAlleleFrequency() >= minimumAlternateAlleleFraction)
						snp.numFilteredHets++;								
				}
			}
			else
				if(debug) System.out.println("empty sample info");
		}
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static HashMap<String, String> getSampleInfoLookup(String formatString, String sampleInfo)
	{
		//format: e.g.  GT:DP:AD:RO:QR:AO:QA
		HashMap<String, String> lookup = new HashMap<String, String>();
		String [] formatTokens = formatString.split(":");
		String [] infoTokens = sampleInfo.split(":");
		for (int i = 0; i < formatTokens.length; i++)
			lookup.put(formatTokens[i], infoTokens[i]);
		
		return lookup;		
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static HashMap<String, String> getInfoLookup(String infoField)
	{
//		System.out.println("infoField = " + infoField);
		HashMap<String, String> lookup = new HashMap<String, String>();
		
		//the info field contains a variable number of key-value pairs separated by semicolons
		String [] infoTokens = infoField.split(";");
		for (int i = 0; i < infoTokens.length; i++)
		{
			//each key-value pair is separated by a "="
			String [] keyValueTokens = infoTokens[i].split("=");
			lookup.put(keyValueTokens[0], keyValueTokens[1]);
		}
		
		return lookup;		
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public LinkedList<VCFEntry> parseVCF(File vcfFile)
	{
		LinkedList<VCFEntry> entries = new LinkedList<VCFEntry>();
		
		try
		{
			int count = 0;
			BufferedReader reader = new BufferedReader(new FileReader(vcfFile));			
			String line = null;
			while((line = reader.readLine()) != null)
			{
				count++;
				if(count % 100000 == 0)
					System.out.print("\rparsing line " + count);
				
				//parse header
				if(line.startsWith("#"))
				{
					//these are the headers that start with two # characters
					if(line.startsWith("##"))
					{
						//add these to a separate list
						doubleHashHeaders.add(line);
					}
					//these are the headers that start with one # character
					else
					{					
						mainHeaders = line.split("\t");
						
						//extract the sample names only
						sampleNamesList = new ArrayList<String>();
						if(debug) System.out.println("extracting headers");
						for (int i = 9; i < mainHeaders.length; i++)
						{
							sampleNamesList.add(mainHeaders[i]);
							if(debug) System.out.println("added header " + mainHeaders[i]);
						}
					}
				}
				//non-comment lines
				else if (!line.startsWith("#") && !line.startsWith("##"))
				{
					VCFEntry vcfEntry = parseVCFLine(line);
					entries.add(vcfEntry);
				}

			}
			
			reader.close();
			
			System.out.println("\ndone parsing VCF file");
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		
		
		
		return entries;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static VCFEntry parseVCFLine(String line)
	{
		//input line looks like this:
		//				#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Barke	Bowman	Morex
		//				TCONS_00000008	225	.	C	A	145	.	DP=13;VDB=0.0518;AF1=1;AC1=6;DP4=0,0,3,10;MQ=20;FQ=-33.4	GT:PL:DP:GQ	1/1:20,3,0:1:12	1/1:20,3,0:1:12	1/1:141,33,0:11:41
		String [] tokens = line.split("\t");
		
		VCFEntry vcfEntry = new VCFEntry();
		vcfEntry.contig = tokens[0].trim();				
		vcfEntry.position = Integer.parseInt(tokens[1]);
		vcfEntry.ID = tokens[2];
		vcfEntry.ref = tokens[3];
		vcfEntry.alt = tokens[4];
		
		//some variant callers, e.g. VarScan, output no QUAL values
		if(!tokens[5].equals("."))
			vcfEntry.qual = Float.parseFloat(tokens[5]);
		else
			vcfEntry.qual = 0;
		
		vcfEntry.filter = tokens[6];
		vcfEntry.info = tokens[7];
		vcfEntry.format = tokens[8];
		
		if(debug) System.out.println("processing SNP on " + vcfEntry.contig + ", pos " + vcfEntry.position );
		
		//check for sample info columns
		if(tokens.length > 8)
		{
			for (int i = 9; i < tokens.length; i++)
			{
				vcfEntry.sampleInfos.add(tokens[i]);
				if(debug) System.out.println("added sample info: " + tokens[i]);
			}
		}	
		
		return vcfEntry;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static ArrayList<String> parseAllVCFHeaders(File vcfFile)
	{
		ArrayList<String> headers = new ArrayList<String>();
		
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(vcfFile));			
			String line = null;
			while((line = reader.readLine()) != null)
			{
				if(line.startsWith("#") || line.startsWith("##"))
					headers.add(line + "\n");
			}
			
			reader.close();
			
			System.out.println("\ndone parsing VCF header");
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		
		return headers;
	}
	
	//---------------------------------------------------------------------------------------------------------------------------------
	
	//This is a cut-down version of parsing a VCF file which is based on just parsing the contig name, the position, 
	//and the alleles representing the variant
	//used for simple tasks that require reading a whole VCF file into memory, e.g. subtraction of variant lists
	public static HashMap<String,VCFEntryLite> parseVCFLite(File vcfFile)
	{
		NumberFormat formatter = new DecimalFormat("###,###");
		
		HashMap<String,VCFEntryLite> entries = new HashMap<String,VCFEntryLite>();
		try
		{
			int count = 0;
			BufferedReader reader = new BufferedReader(new FileReader(vcfFile));			
			String line = null;
			while((line = reader.readLine()) != null)
			{
				count++;
//				if(count % 100 == 0)
//					System.out.print("\rparsing line " + formatter.format(count));
				
				//only the non-commented lines contain the actual entries
				//we are not interested in headers for this purpose
				if (!line.startsWith("#") && !line.startsWith("##"))
				{
					VCFEntryLite vcfEntry = parseVCFLineLite(line);
					//now store the entry in the lookup, using a concatenation of the contig name and the position as key, and
					//the object itself as the value
					String key = vcfEntry.contig+"_"+vcfEntry.position;
					entries.put(key,vcfEntry);
				}
			}
			
			reader.close();
			
			System.out.println("done parsing VCF file");
			System.out.println("# variants parsed = " + formatter.format(count));
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
			
		return entries;
	}
	
	//------------------------------------------------------------------------------------------------------------------------

	//This is a cut-down version of parsing a VCF file which is based on just parsing the contig name, the position, 
	//and the alleles representing the variant
	//used for simple tasks that require reading a whole VCF file into memory, e.g. subtraction of variant lists
	public static VCFEntryLite parseVCFLineLite(String line)
	{
		//input line looks like this:
		//				#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Barke	Bowman	Morex
		//				TCONS_00000008	225	.	C	A	145	.	DP=13;VDB=0.0518;AF1=1;AC1=6;DP4=0,0,3,10;MQ=20;FQ=-33.4	GT:PL:DP:GQ	1/1:20,3,0:1:12	1/1:20,3,0:1:12	1/1:141,33,0:11:41
		String [] tokens = line.split("\t");
		
		VCFEntryLite vcfEntry = new VCFEntryLite();
		vcfEntry.contig = tokens[0].trim();				
		vcfEntry.position = Integer.parseInt(tokens[1]);
		vcfEntry.ID = tokens[2];
		vcfEntry.ref = tokens[3];
		vcfEntry.alt = tokens[4];
		
		//also add a corresponding SNP object
		vcfEntry.snp = new SNP();
		vcfEntry.snp.refSeqName = vcfEntry.contig;
		vcfEntry.snp.position = vcfEntry.position;
		vcfEntry.snp.refSeqAllele = new Allele(vcfEntry.ref.charAt(0));
		vcfEntry.snp.alternateAllele = new Allele(vcfEntry.alt.charAt(0));
		
		
		if(debug) System.out.println("processing SNP on " + vcfEntry.contig + ", pos " + vcfEntry.position );
	
		
		return vcfEntry;
	}
}
