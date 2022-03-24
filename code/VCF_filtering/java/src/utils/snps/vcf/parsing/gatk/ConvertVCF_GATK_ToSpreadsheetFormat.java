package utils.snps.vcf.parsing.gatk;

import java.io.*;
import java.text.*;
import java.util.*;
import utils.entities.*;
import utils.fasta.*;
import utils.snps.SNPUtils;
import utils.snps.filtering.filters.*;
import utils.snps.vcf.*;
import utils.snps.vcf.parsing.*;
import utils.validation.transcripts.*;


/**
 * Takes VCF output from GATK and parses it. Outputs a human-friendly spreadsheet format that can be used to 
 * readily filter SNPs as part of the downstream processing.   
 */
public class ConvertVCF_GATK_ToSpreadsheetFormat
{
	static boolean debug = false;	

	final static String USAGE = "java utils.snps.vcf.parsing.gatk.ConvertVCF_GATK_ToSpreadsheetFormat <VCF file> <config file>" ;
	
	public static File vcfFile;
	
	public static String [] headers;
	public static ArrayList<String> sampleNamesList;	

	//parameters for filtering
	public static float minorGenotypeFrequencyCutoff;
	public static int minorGenotypeCountCutoff;
	public static int qualScoreCutoff;
	public static float pctHetSamplesCutoff;
	public static int minCoverageDepth;
	public static float minCoverageSamplesPercentage;
	public static int ploidy =2;
	public static float missingSamplesCutoff;

	//a list of the corresponding filter objects
	static ArrayList<SNPFilter> filters;
	
	//file writers
	static BufferedWriter filteredSNPsWriter,  filteredGffWriter,  rejectedSNPsWriter, rejectedGffWriter, snpVCFWriter, filteredIndelsVCFWriter, filteredIndelsTxtWriter, mnpVCFWriter;
	
	//formatter for the count output
	static NumberFormat formatter = new DecimalFormat("###,###");

	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------	
	
	public static void main(String[] args)
	{
		if(args.length != 2)
		{
			System.out.println("\nERROR: incorrect number of args supplied.\nCorrect usage is:" + USAGE);
			System.exit(1);
		}
		
		try
		{
			System.out.println("ConvertVCF_GATK_ToSpreadsheetFormat");
			System.out.println("starting filtering run");
			
			ConvertVCF_GATK_ToSpreadsheetFormat converter = new ConvertVCF_GATK_ToSpreadsheetFormat();
			
			//the vcfFile we need to process
			vcfFile = new File(args[0]);
			//its base name
			String baseName = vcfFile.getName().substring(0, vcfFile.getName().indexOf("."));
			System.out.println("VCF file for filtering = " + vcfFile.getName());
			
			//parse the config file and configure all required parameters
			configureParameters(new File(args[1]), baseName);
			
			//make a list of filter objects
			//these same filters will also be added to each SNP but for convenience we also hold the same list here
			//a list of the corresponding filter objects
			filters = addFilters();

			//go process
			processFile();
			
			//close file writers
			converter.closeWriters();
			
			System.out.println("\ndone");
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private void closeWriters() throws IOException
	{
		rejectedSNPsWriter.close();
		filteredSNPsWriter.close();
		filteredGffWriter.close();
		rejectedGffWriter.close();
		snpVCFWriter.close();
		filteredIndelsVCFWriter.close();
		mnpVCFWriter.close();
		filteredIndelsTxtWriter.close();
	}
	

	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static ArrayList<SNPFilter> addFilters()
	{
		//filtering required:
		
//		MIN_COVERAGE=8,50
//		MISSING_SAMPLES_CUTOFF=5
//		MINOR_GENOTYPE_FREQ=5
//		MINOR_GENOTYPE_COUNT=1
//		QUAL=30
//		PCT_HET_SAMPLES_CUTOFF=2
//		INDEL_REMOVAL=TRUE
				
		//a list of the corresponding filter objects
		ArrayList<SNPFilter> filters = new ArrayList<SNPFilter>();
		
		//add the appropriate filters we want and configure each of them with the appropriate parameters
		filters.add(new QUALFilter(qualScoreCutoff));
		filters.add(new MinorGenotypeFreqFilter(minorGenotypeFrequencyCutoff));
		filters.add(new MinorGenotypeCountFilter(minorGenotypeCountCutoff));
		filters.add(new HeterozygousSampleFilter(pctHetSamplesCutoff));
		filters.add(new SampleDepthFilter(minCoverageDepth,minCoverageSamplesPercentage));
		filters.add(new MissingDataFilter(missingSamplesCutoff));
		
		return filters;
	}

	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static HashMap<String,String> parseConfigFile(File configFile) throws IOException
	{
		System.out.println("parsing config file");
		
		HashMap<String,String> parameterMap = new HashMap<String,String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(configFile));			
		String line = null;

		//parse the file
		while((line = reader.readLine()) != null)
		{
			//ignore header lines 
			if (!line.startsWith("#") && !line.startsWith(" ") && line!=null && !line.equals(""))
			{
				if(debug) System.out.println("config file line: " + line);
				
				String [] tokens = line.split("=");
				parameterMap.put(tokens[0], tokens[1]);
			}
		}
		
		reader.close();
		
		return parameterMap;
	}
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static void configureParameters(File configFile, String baseName) throws IOException 
	{
		//read the config file and return a String based map of params from it
		HashMap<String,String> parameterMap = parseConfigFile(configFile);

		System.out.println("list of parameter values supplied in config file:");
		for(String param: parameterMap.keySet())
			System.out.println(param + "\t" + parameterMap.get(param));
		
		//sample config file content and format:
//		VCF_PATH=
//		MIN_COVERAGE=8,50
//		MISSING_SAMPLES_CUTOFF=5
//		MINOR_GENOTYPE_FREQ=5
//		QUAL=30
//		PCT_HET_SAMPLES_CUTOFF=2
		
//		the minum coverage required to retain a variant
//		this is a composite value
//		specified as e.g. "8,50" which means we want at least 8x coverage for at least 50% of the samples
		String minCoverage = parameterMap.get("MIN_COVERAGE");
		String [] covTokens = minCoverage.split(",");
		minCoverageDepth = Integer.parseInt(covTokens[0]);
		minCoverageSamplesPercentage = Float.parseFloat(covTokens[1]);
		
//		a threshold for discarding a variant based on samples having no data
//		e.g. a value of "5" indicates that if more than 5% of samples have no data the variant is removed
		missingSamplesCutoff = Float.parseFloat(parameterMap.get("MISSING_SAMPLES_CUTOFF"));

		//the minor genotype frequency in the population of samples (not among reads)
		//e.g. a value of 5 means we want at least 5% of samples to have the minor genotype to retain the variant
		minorGenotypeFrequencyCutoff = Float.parseFloat(parameterMap.get("MINOR_GENOTYPE_FREQ"));
		
		//the minimum number of samples required with the minor genotype
		//samples have to be homozygous for this
		//this discards SNPs that have been triggered by one or more heterozygous samples only
		minorGenotypeCountCutoff = Integer.parseInt(parameterMap.get("MINOR_GENOTYPE_COUNT"));	
		
		//the min. quality score to retain a variant
		qualScoreCutoff = Integer.parseInt(parameterMap.get("QUAL"));
		
		//the percentage cutoff of samples that have been called as heterozygotes
		//if the variant has a percentage higher than this threshold, it will be discarded
		//this filter only works for species where we are expecting very little heterozygosity and whatever little there is, is expected to have different genomic distribution patterns based in every sample
		//this means we can remove SNPs that are caused by mismapping, if we suddenyl have lots of samples that are heterozygous at the same location, which is biologically implausible
		pctHetSamplesCutoff = Float.parseFloat(parameterMap.get("PCT_HET_SAMPLES_CUTOFF"));
		
		//output files
		
		//GFF file with the good (filtered) SNPs
		File filteredGffFile = new File(baseName + "_filteredSNPs.gff");
		filteredGffWriter = new BufferedWriter(new FileWriter(filteredGffFile));
		
		//GFF file  with the bad (rejected) SNPs
		File rejectedGffFile = new File(baseName + "_rejectedSNPs.gff");
		rejectedGffWriter = new BufferedWriter(new FileWriter(rejectedGffFile));
		
		//file with the good (filtered) SNPs
		File filteredSNPsFile = new File(baseName + "_filteredSNPs.txt");
		filteredSNPsWriter = new BufferedWriter(new FileWriter(filteredSNPsFile));
		
		//file with the bad (rejected) SNPs
		File rejectedSNPsFile = new File(baseName + "_rejectedSNPs.txt");
		rejectedSNPsWriter = new BufferedWriter(new FileWriter(rejectedSNPsFile));	
		
		//VCF file with the filtered (accepted) SNPs
		File vcfOutputFile = new File(baseName + "_filteredSNPs.vcf");
		snpVCFWriter = new BufferedWriter(new FileWriter(vcfOutputFile));	
		
		//VCF and .txt  files with the filtered (accepted) indels
		File indelVCFFile = new File(baseName + "_filteredIndels.vcf");
		File indelsTxtFile =  new File(baseName + "_filteredIndels.txt");
		filteredIndelsTxtWriter = new BufferedWriter(new FileWriter(indelsTxtFile));				
		filteredIndelsVCFWriter = new BufferedWriter(new FileWriter(indelVCFFile));	
		
		//VCF file with the rejected MNPs and multiallelic SNPs
		File mnpVCFFile = new File(baseName + "_mnps.vcf");
		mnpVCFWriter = new BufferedWriter(new FileWriter(mnpVCFFile));	
	}

	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static void processFile() throws IOException
	{
		//parse the VCF header 
		ArrayList<String> headerList = parseHeader(vcfFile);
		
		//then write headers for the output files
		SNPOutput.writeHeaderToFile(filteredSNPsWriter, sampleNamesList);
		SNPOutput.writeHeaderToFile(rejectedSNPsWriter, sampleNamesList);
		SNPOutput.writeHeaderToFile(filteredIndelsTxtWriter, sampleNamesList);
		VCFEmitter.outputVCFHeader(snpVCFWriter, headerList);
		VCFEmitter.outputVCFHeader(filteredIndelsVCFWriter, headerList);
		VCFEmitter.outputVCFHeader(mnpVCFWriter, headerList);
		
		System.out.println("processing VCF file");
		
		//I/O
		BufferedReader reader = new BufferedReader(new FileReader(vcfFile));			
		String line = null;
		
		//counts
		int count = 0;
		int numSNPsFailedFilter = 0;
		int numSNPsOutput = 0;
		int numMultiAllelicSNPs = 0;
		int numIndelsOutput = 0;
		
		//parse the file
		while((line = reader.readLine()) != null)
		{
			//ignore header lines 
			if (!line.startsWith("#") && !line.startsWith("##"))
			{
				count++;
				if(count % 1000 == 0)
					System.out.print("\rprocessing SNP " + formatter.format(count));
						
				//make a new SNP object
				SNP snp = new SNP();
				//set the line as a variable on this
				snp.vcfLine = line;		
				//configure the new SNP object appropriately
				boolean success = SNPConfigurator.configureSNP(snp, debug, sampleNamesList, ploidy, minCoverageDepth);	
				
				//check the SNP was configured ok
				//if it wasn't, it is either an MNP ("complex" type) or multiallelic SNP
				//these get discarded, but we still write them to a (separate) VCF file for checking
				if(!success)
				{
					SNPOutput.outputMNP(snp);
					numMultiAllelicSNPs++;
					continue;
				}
				
				//check the SNP against the filtering criteria passed in as parameters 
				boolean pass = SNPConfigurator.filterSNP(snp, debug);	
				//this is a fail
				if(!pass)
				{
					numSNPsFailedFilter++;	
				}
				//this is a pass
				else
				{
					//for proper SNPs, just increment the count here
					if(!snp.isIndel)
						numSNPsOutput++;
					//count and output indels separately, and only those that passed
					else
					{		
						SNPOutput.outputIndel(snp);
						numIndelsOutput++;
					}
				}

				//output all SNPs to file regardless of whether they fail or pass, but separate passes and fails into separate files for checking
				if(!snp.isIndel)
					SNPOutput.outputSNP(snp, pass);	
						
			}			
		}
		
		//done -- close the reader
		reader.close();
		
		//now print out some stats
		System.out.println();
		System.out.println("# variants processed: " + formatter.format(count));
		System.out.println("# SNPs passed and output to file: " + formatter.format(numSNPsOutput));
		System.out.println("# indels passed and output to file: " + formatter.format(numIndelsOutput));
		System.out.println("# MNPs and multiallelic variants discarded: " + formatter.format(numMultiAllelicSNPs));
		System.out.println("# SNPs/indels discarded due to failed filters: " + formatter.format(numSNPsFailedFilter));
		float pctDiscarded = numSNPsFailedFilter/(float)count;		
		NumberFormat percentFormat = NumberFormat.getPercentInstance();
		percentFormat.setMaximumFractionDigits(2);
		System.out.println("% SNPs/indels discarded due to failed filters: " + percentFormat.format(pctDiscarded));
		System.out.println("\ndone parsing VCF file");

	}

	
	//-----------------------------------------------------------------------------------------------------------------------------------------------
	
	public static ArrayList<String> parseHeader(File vcfFile) throws IOException
	{
		System.out.println("parsing header");

		ArrayList<String> headerList = new ArrayList<String>();
			
		BufferedReader reader = new BufferedReader(new FileReader(vcfFile));			
		String line = null;
		
		boolean headersParsed = false;
		while((line = reader.readLine()) != null && !headersParsed)
		{
			//parse header
			if(line.startsWith("#"))//all headers start with at least one #
			{
				if(!line.startsWith("##"))//the main header only has a single # prefix
				{
					headerList.add(line+"\n");
					
					headers = line.split("\t");
					
					//extract the sample names only
					sampleNamesList = new ArrayList<String>();
					if(debug) System.out.println("extracting headers");
					for (int i = 9; i < headers.length; i++)
					{
						sampleNamesList.add(headers[i]);
						if(debug) System.out.println("added header " + headers[i]);
					}	
					
					headersParsed = true;
				}
				else//all the other headers have two # as a prefix
				{
					headerList.add(line+"\n");
				}
			}
		}
		reader.close();
		
		System.out.println("done parsing header");
		
		return headerList;
		
	}

	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

}
