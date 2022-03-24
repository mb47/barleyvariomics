package utils.snps.vcf.parsing.gatk;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import utils.entities.*;
import utils.snps.filtering.filters.SNPFilter;
import utils.snps.vcf.parsing.VCFEmitter;


public class SNPOutput
{
	
	static DecimalFormat percentageFormatter = new DecimalFormat("#.##");	

	//------------------------------------------------------------------------------------------------------------------------		
	
	public static void outputSNP(SNP snp, boolean passedFilters) throws IOException
	{
		//if it has passed, add it to the list of SNPs returned
		if(passedFilters)
		{
			//output SNP
			SNPOutput.writeSNPInfoToFile(ConvertVCF_GATK_ToSpreadsheetFormat.filteredSNPsWriter, snp, ConvertVCF_GATK_ToSpreadsheetFormat.sampleNamesList);
			SNPOutput.writeSNPInfoToGFF(ConvertVCF_GATK_ToSpreadsheetFormat.filteredGffWriter, snp, "GATK", "SNP_ACCEPTED");
			VCFEmitter.outputSingleVCFEntry(ConvertVCF_GATK_ToSpreadsheetFormat.snpVCFWriter, snp.vcfEntry);
		}
		//else output it to the file with the rejected SNPs
		else
		{
			SNPOutput.writeSNPInfoToFile(ConvertVCF_GATK_ToSpreadsheetFormat.rejectedSNPsWriter, snp, ConvertVCF_GATK_ToSpreadsheetFormat.sampleNamesList);
			SNPOutput.writeSNPInfoToGFF(ConvertVCF_GATK_ToSpreadsheetFormat.rejectedGffWriter, snp, "GATK", "SNP_REJECTED");				
		}
	}
	
	//------------------------------------------------------------------------------------------------------------------------		
	
	public static void outputMNP(SNP snp) throws IOException
	{
		//output indel -- we only get this far if it passed all filters anyway
		VCFEmitter.outputSingleVCFEntry(ConvertVCF_GATK_ToSpreadsheetFormat.mnpVCFWriter, snp.vcfEntry);
	}
	
	//------------------------------------------------------------------------------------------------------------------------		
	
	public static void outputIndel(SNP indel) throws IOException
	{
		//output indel -- we only get this far if it passed all filters anyway
		VCFEmitter.outputSingleVCFEntry(ConvertVCF_GATK_ToSpreadsheetFormat.filteredIndelsVCFWriter, indel.vcfEntry);
		
		SNPOutput.writeSNPInfoToFile(ConvertVCF_GATK_ToSpreadsheetFormat.filteredIndelsTxtWriter, indel, ConvertVCF_GATK_ToSpreadsheetFormat.sampleNamesList);
	}
	
	//------------------------------------------------------------------------------------------------------------------------		

	public static void writeHeaderToFile(BufferedWriter writer, List<String> sampleNamesList) throws IOException
	{
		writer.write("SNP name" + "\t");
		writer.write("reference sequence name" + "\t");
		writer.write("SNP position" + "\t");
		
		writer.write("quality score" + "\t");
	
		writer.write("reference allele" + "\t");
		writer.write("alternate allele" + "\t");
		writer.write("major allele" + "\t");
		writer.write("minor allele" + "\t");

		writer.write("total read coverage" + "\t");	
		writer.write("# samples in dataset" + "\t");
		writer.write("# samples with data in SNP" + "\t");
		
		writer.write("# heterozygous samples in SNP" + "\t");
		writer.write("% heterozygous samples in SNP" + "\t");
	
		writer.write("# samples with major allele as genotype" + "\t");
		writer.write("# samples with minor allele as genotype" + "\t");
		writer.write("% samples with major allele as genotype" + "\t");
		writer.write("% samples with minor allele as genotype" + "\t");
		
		writer.write("% samples with required coverage (" + ConvertVCF_GATK_ToSpreadsheetFormat.minCoverageDepth + "x)\t");
		
		writer.write("% samples with no data" + "\t");

		// list all filters; pass/fail (true/false)
		for(SNPFilter filter: ConvertVCF_GATK_ToSpreadsheetFormat.filters)
			writer.write(filter.label + "\t");
		
		for (String sampleName : sampleNamesList)
		{
			writer.write(sampleName + "\t");
		}
	
		writer.write("read counts all samples" + "\t");
		
		writer.write("SnpEff annotation" + "\t");
		
		writer.newLine();
	}
	
	//------------------------------------------------------------------------------------------------------------------------	

	public static void writeSNPInfoToGFF(BufferedWriter writer, SNP snp, String source, String type) throws IOException
	{
		writer.write(snp.refSeqName + "\t"+source+"\t" + type + "\t");
		writer.write(snp.position + "\t" + snp.position + "\t.\t.\t.\tName=" + snp.name + "\n");
	}
	
	//------------------------------------------------------------------------------------------------------------------------	

	public static void writeSNPInfoToFile(BufferedWriter writer, SNP snp, List<String> sampleNamesList) throws IOException
	{
		String outputLine = prepareLineForOutput(snp,sampleNamesList);
		writer.write(outputLine);
		writer.newLine();
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	public static String prepareLineForOutput(SNP snp, List<String> sampleNamesList)
	{
		StringBuilder builder = new StringBuilder(); 
		
		builder.append(snp.name + "\t");
		builder.append(snp.refSeqName + "\t");
		builder.append(snp.position + "\t");

		builder.append(snp.qualityScore + "\t");
		
		builder.append(snp.refSeqAllele + "\t");
		builder.append(snp.alternateAllele + "\t");
		builder.append(snp.majorAllele + "\t");
		builder.append(snp.minorAllele + "\t");
			
		builder.append(snp.totalReadCount + "\t");
		builder.append(sampleNamesList.size() + "\t");
		builder.append(snp.samplesList.size() + "\t");
		
		builder.append(snp.numHets + "\t");		
		// same but percentages
		if (snp.samplesList.size() > 0)
		{
			builder.append(percentageFormatter.format(snp.pctHets) + "\t");
		}
		else
		{
			builder.append("0\t");
		}
		
		// number of samples that have the major/minor allele as their genotype
		builder.append(snp.numberOfSamplesWithMajorAlleleAsGenotype + "\t");
		builder.append(snp.numberOfSamplesWithMinorAlleleAsGenotype + "\t");
//		the same as percentages
		builder.append(percentageFormatter.format(snp.pctSamplesWithMajorAlleleAsGenotype) + "\t");
		builder.append(percentageFormatter.format(snp.pctSamplesWithMinorAlleleAsGenotype) + "\t");

		//percentage of samples that have at least the level of coverage required by our filter
		builder.append(snp.pctSamplesWithSufficientDepth + "\t");
		
		//percentage of samples with no data
		builder.append(snp.pctMissingSamples + "\t");
		
		// list all filters; pass/fail (true/false)
		for(SNPFilter filter: snp.filters)
		{
			//check how we need to format the output value for this filter
			if(filter.formatOutputValueAsPercent)
				builder.append(filter.getPassMarkAsString() + " (" + percentageFormatter.format(filter.getTestedValue()) + ")" + "\t");
			else
				builder.append(filter.getPassMarkAsString() + " (" + filter.getTestedValue() + ")" + "\t");
		}
		
		// sample genotypes
		for (String sampleName : sampleNamesList)
		{
			SNPSample sample = snp.getSNPSampleByName(sampleName);
			if (sample != null)
			{
				builder.append(sample.genotypeFromVCFAsString + "\t");
			}
			else
				builder.append("-\t");
		}
		
		// this outputs allele counts for all samples in a single column, separated by "|"
		for (String sampleName : sampleNamesList)
		{
			SNPSample snpSample = snp.getSNPSampleByName(sampleName);		
			if (snpSample != null)
			{
				builder.append(sampleName + " : ");
				
				for (Allele allele : snpSample.alleleList)
				{
					builder.append(" " + allele.stringValue + " " + allele.count + ";");
				}
				builder.append("|");
			}
		}
		
		//output the SNP effect annotation string if provided by snpEff
		if(snp.snpEffAnotationString != null)
			builder.append("\t" + snp.snpEffAnotationString + "\t");

		return builder.toString();
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
}
