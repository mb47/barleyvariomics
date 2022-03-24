package utils.snps.vcf.parsing.gatk;

import java.util.*;
import utils.entities.*;
import utils.snps.SNPUtils;
import utils.snps.filtering.filters.SNPFilter;
import utils.snps.vcf.parsing.VCFParser;

public class SNPConfigurator
{

	
	//------------------------------------------------------------------------------------------------------------------------	

	public static boolean configureSNP(SNP snp, boolean debug, ArrayList<String> sampleNamesList, int ploidy, int readDepthCutoff)
	{
		//parse the VCF line into in a VCFEntry object
		snp.vcfEntry = VCFParser.parseVCFLine(snp.vcfLine);		
		
		//configure this SNP with the info from the VCF entry
		snp.refSeqName = snp.vcfEntry.contig;
		snp.position = snp.vcfEntry.position;
		snp.name = snp.refSeqName + "_" + snp.position;
		snp.qualityScore = snp.vcfEntry.qual;
		
		if(debug) System.out.println("\n=========parsing SNP at pos " + snp.position + " in contig " + snp.refSeqName);
		
		//deal with MNPs ("complex"), indels and multiallelic SNPs
		
		//multi-allelic SNPs
		//example: alt is G,C
		snp.isMultiAllelic = snp.vcfEntry.alt.contains(",");
		
		//indels
		//example: ref is A, alt is AGAG (but alt can't contain multiple alleles) 
		snp.isIndel = snp.vcfEntry.ref.length() > 1 || (snp.vcfEntry.alt.length() > 1 && !snp.isMultiAllelic);
		
		//multi-nucleotide polymorphisms
		//example: ref is AGA, alt is GAA
		snp.isMNP = (snp.vcfEntry.ref.length() == snp.vcfEntry.alt.length()) && !snp.isMultiAllelic && snp.vcfEntry.ref.length() > 1;
		
		if(debug) 
		{
			System.out.println("snp.vcfEntry.ref = " + snp.vcfEntry.ref);
			System.out.println("snp.vcfEntry.alt = " + snp.vcfEntry.alt);
			System.out.println("snp.isMultiAllelic = " + snp.isMultiAllelic);
			System.out.println("snp.isIndel = " + snp.isIndel);
			System.out.println("snp.isMNP = " + snp.isMNP);
		}
		
		//we throw out anything that is multi-allelic or MNP
		if(snp.isMultiAllelic || snp.isMNP)
		{	
			return false;
		}
		
		//assign alleles
		snp.refSeqAllele = new Allele(snp.vcfEntry.ref);
		snp.alternateAllele = new Allele(snp.vcfEntry.alt);
		
		if(debug)
		{
			System.out.println("snp.refSeqAllele = " + snp.refSeqAllele);
			System.out.println("snp.alternateAllele = " + snp.alternateAllele);
		}
		
		//parse info string and put into map
		HashMap<String, String> infoMap = SNPUtils.parseInfoString(snp.vcfEntry);
		
		//check for SNP effect information provided by snpEff
		//this uses  the "EFF" tag
		String snpEffInfo = infoMap.get("EFF");
		if(snpEffInfo != null)
			snp.snpEffAnotationString = snpEffInfo;

		snp.totalReadCount = Integer.parseInt(infoMap.get("DP"));
		
		//populate list of samples for this SNP
		if(debug) System.out.println("vcfEntry.sampleInfos.size() = " + snp.vcfEntry.sampleInfos.size());		
		processSampleInfos(snp.vcfEntry, snp, sampleNamesList, debug);	
		
		//compute the number of samples that are homozygous for either the ref or alt allele
		//number of samples that have the major/minor allele as their genotype
//		System.out.println("\n\ncounting samples with ref genotype");
		int numSamplesWithRefGT = snp.computeNumberOfSamplesWithAlleleAsGenotype(snp.refSeqAllele, ploidy);
//		System.out.println("\n\ncounting samples with alt genotype");
		int numSamplesWithAltGT = snp.computeNumberOfSamplesWithAlleleAsGenotype(snp.alternateAllele, ploidy);
				
		//assign the minor and major allele status accordingly
		if(numSamplesWithRefGT >= numSamplesWithAltGT)	
		{
			snp.majorAllele = snp.refSeqAllele;
			snp.minorAllele = snp.alternateAllele;
			snp.numberOfSamplesWithMajorAlleleAsGenotype = numSamplesWithRefGT;
			snp.numberOfSamplesWithMinorAlleleAsGenotype = numSamplesWithAltGT;
		}
		else
		{
			snp.majorAllele = snp.alternateAllele;
			snp.minorAllele = snp.refSeqAllele;
			snp.numberOfSamplesWithMajorAlleleAsGenotype = numSamplesWithAltGT;
			snp.numberOfSamplesWithMinorAlleleAsGenotype = numSamplesWithRefGT;
		}
		
		if(debug)
		{		
			System.out.println("numSamplesWithRefGT = " + numSamplesWithRefGT);
			System.out.println("numSamplesWithAltGT = " + numSamplesWithAltGT);
			System.out.println("snp.majorAllele = " + snp.majorAllele);
			System.out.println("snp.minorAllele = " + snp.minorAllele);
		}
		
		//work out what this equates to as percentages
		snp.pctSamplesWithMajorAlleleAsGenotype = (snp.numberOfSamplesWithMajorAlleleAsGenotype / (float)snp.samplesList.size()) * 100;
		snp.pctSamplesWithMinorAlleleAsGenotype = (snp.numberOfSamplesWithMinorAlleleAsGenotype / (float)snp.samplesList.size()) * 100;
		
		//work out % of heterozygous samples
		snp.pctHets = ((snp.numHets/(float)snp.samplesList.size())*100);
		
		//now add a list of filter objects to this, in preparation of the filtering later
		snp.filters = ConvertVCF_GATK_ToSpreadsheetFormat.addFilters();
		
		//this is the total number of samples in our dataset, i.e. the max. achievable
		int totalNumSamplesInDataset = sampleNamesList.size();
		//this is the number of samples that actually have ANY data in this current SNP
		int numSamplesInSNP = snp.samplesList.size();
		//work out how much is missing
		int numMissingSamples = totalNumSamplesInDataset - numSamplesInSNP;
		//and the same as a percentage
		snp.pctMissingSamples = (numMissingSamples / (float) totalNumSamplesInDataset) *100;
		
		//work out the number of samples that have at least the read depth we want to filter on
		//the total number of samples
		int numSamples = snp.samplesList.size();		
		//a counter for the number of samples that exceed our required threshold
		int numSamplesWithSufficientDepth = 0;		
		//iterate over the samples and chek them individually
		for(SNPSample sample: snp.samplesList)
		{
			if(sample.readCount >= readDepthCutoff)
				numSamplesWithSufficientDepth++;
		}		
		//now work out the percentage of samples that satisfy this requirement
		snp.pctSamplesWithSufficientDepth = (numSamplesWithSufficientDepth / (float)numSamples) * 100;
		
		
		return true;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	static void processSampleInfos(VCFEntry vcfEntry, SNP snp, ArrayList<String> sampleNamesList, boolean debug)
	{
		//format of sample info:
		/*
		GT:AD:DP:GQ:PGT:PID:PL
		0/0:3,0:3:9:.:.:0,9,79
		
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
		##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
		##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,D	escription="Minimum DP observed within the GVCF block">
		##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
		##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
		##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
		 */
		
		//populate list of samples for this SNP
		//iterate over all sample infos
		//if the text is not "./.:0,0:0" then we add the sample to our list
		for(int i =0; i < vcfEntry.sampleInfos.size(); i++)
		{
			String sampleInfo = vcfEntry.sampleInfos.get(i);
			String sampleName =  sampleNamesList.get(i);
			
			if(debug) System.out.println("\n+++++sampleIndex = " + i + ", sample " + sampleName);
			
			if(!sampleInfo.startsWith("."))
			{
				//get the name of this sample
				SNPSample snpSample = new SNPSample();
				snpSample.name = sampleName;

				//extract sample info
				HashMap<String, String> lookup = SNPUtils.getSampleInfoLookup(vcfEntry.format, sampleInfo);
				
				//extract the read depth for this sample
				String dpString = lookup.get("DP");
				//we need to check that the DP string is not "." - this means the read depth for this sample is zero, but
				//for some reason GATK has still assigned a genotype to it (presumably based on interpolation?)
				//in this case we don't want to retain this sample
				if(!dpString.equals("."))
				{
					//add the sample to the list
					snp.samplesList.add(snpSample);
					snpSample.readCount = Integer.parseInt(dpString);
				}
				//skip to the next sample if this one has a DP of 0
				else
					continue;
				
				//extract allele counts and add to list 
				String adString = lookup.get("AD");
				String [] adTokens = adString.split(",");
				//check if we have more than one alt allele -- this shouldn't happen
				if(adTokens.length > 2)
				{
					System.out.println("ERROR: more than one alt allele found -- exiting.");
					System.exit(1);
				}
				//also add the alleles to the sample's allele list for convenience
				snpSample.refSeqAllele = new Allele (snp.refSeqAllele.stringValue);
				snpSample.alternateAllele = new Allele (snp.alternateAllele.stringValue);			
				snpSample.alleleList.add(snpSample.refSeqAllele);
				snpSample.alleleList.add(snpSample.alternateAllele);
				
				//assign the counts
				snpSample.refSeqAllele.count = Integer.parseInt(adTokens[0]);
				snpSample.alternateAllele.count =  Integer.parseInt(adTokens[1]);	
				
				//extract and convert the VCF genotype
				snpSample.numericalGenotypeFromVCFString = lookup.get("GT");
				snpSample.singleDigitNumericalGenotypeFromVCFString = SNPUtils.convertNumericalVCFGenotypeToSingleDigit(snpSample.numericalGenotypeFromVCFString);
				snpSample.genotypeFromVCFAsString = SNPUtils.translateVCFGenotype(snpSample.numericalGenotypeFromVCFString, snp);	
				//and extract the phred score for this genotype
				//however, this can be null if we run freebayes with the --no-population-priors flag
				//this suppresses the computation of a genotype quality score
				String gtScoreStr = lookup.get("GQ");
				if(gtScoreStr != null)
					snpSample.genotypeFromVCFPhredScore = Float.parseFloat(gtScoreStr);	

				if (debug)
				{
					System.out.println("snpSample.isHomozygous() = " + snpSample.vcfGenotypeIsHomozygous());
					System.out.println("sampleInfo = " + sampleInfo);
				}
				
				//ditto for the smallest USABLE allele count in this sample
				//that is, for homozygotes we use the major count, for hets the minor count
				//this also counts the number of hets and homs
				if (snpSample.vcfGenotypeIsHomozygous())
					snp.numCleanHomoZygotes++;
				else
					snp.numHets++;						
			}
			else
				if(debug) System.out.println("empty sample info");
		}
	}


	//------------------------------------------------------------------------------------------------------------------------

		//runs a SNP through a series of filters and returns a boolean that indicates whether it passed
		public static boolean filterSNP(SNP snp, boolean debug)
		{
			if(debug) System.out.println("\nfiltering SNP " + snp.name);
			
			//first we need to apply all filters
			for(SNPFilter filter: snp.filters)
			{
				filter.filterSNP(snp);
				if(debug) System.out.println("filter " + filter.label + " = " + filter.getPassMarkAsString());
			}
			
			//now iterate over all filters again
			//if even just one of them fails, fail the SNP as a whole
			boolean pass = true;
			for(SNPFilter filter: snp.filters)
			{
				if(!filter.pass)
					pass = false;
			}
			
			if(debug) System.out.println("overall pass for this SNP = " + pass);
			
			return pass;
		}
		
		//------------------------------------------------------------------------------------------------------------------------
}
