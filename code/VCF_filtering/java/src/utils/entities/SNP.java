package utils.entities;

import java.util.*;
import utils.entities.*;
import utils.snps.SNPUtils;
import utils.snps.filtering.filters.SNPFilter;

public class SNP
{
	public String name;
	public String description;
	public String refSeqName;
	public int position;
	public String sequence;
	public SNPSequence snpSequence;
	
	//the corresponding entry in a VCF file
	public VCFEntry vcfEntry;
	//the VCF row of data for this SNP - verbatim as in the original VCF file, no extra punctuation or parsing etc (single line of text)
	public String vcfLine;	
	
	public String comment;
	
	public String lhFlankingRegion;
	public String rhFlankingRegion;
	
	public int numCleanHomoZygotes = 0;
	public int numHets = 0;
	//this is a value we use for getting a better idea of genuine hets
	//for this we filter the allele counts within a given SNP x sample combination the same way we filter the overall allele frequencies
	//i.e. if we filter SNP below a certain minor allele frequency and minro allele count then we do the same within samples here
	public int numFilteredHets = 0;
	//the same as percentages
	public float pctHets = -1;
	public float pctFilteredHets = -1;
	
	public ArrayList<Character> alleles = new ArrayList<Character>();
	public ArrayList<Integer> alleleCounts;
	
	//this is a list of the actual samples that contribute to this SNP (can be less than the total number of samples involved in this run)
	public ArrayList<SNPSample> samplesList = new ArrayList<SNPSample>();
	
	//an equivalent lookup data structure where we store a SNPSample against its name
	public HashMap<String, SNPSample> samplesLookup = new HashMap<String, SNPSample>();
	
	//this lookup holds the number of alleles against the allele char
	public HashMap<Character, Integer> alleleCountsLookup = new HashMap<Character, Integer>();
	
	public int distanceToEndOfReference = -1;
	public boolean refSeqExcluded = false;
	public boolean otherSNPsInFlankingRegion = false;
	public boolean nextSNPIsInFlankingRegion = false;
	
	public ArrayList<SNP> snpsInLeftFlankingRegion = new ArrayList<SNP>();
	public ArrayList<SNP> snpsInRightFlankingRegion = new ArrayList<SNP>();
	
	public int refSeqLength = -1;
	
	//this helps us keep tab of the total number of reads that contribute to the SNP
	public int totalReadCount = 0;
	
	//the minimum total number of reads we have found in any sample for this SNP
	public int minTotalAlleleCountForASample = Integer.MAX_VALUE;
	
	//the smallest number of reads we have found for any allele in any sample for this SNP
	public int minAlleleCountForASample = Integer.MAX_VALUE;
	
	//the smallest ratio of minor allele to total allele count we have found in any sample for this SNP
	public float minAlleleRatioForASample = Float.MAX_VALUE;
	
	public Allele minorAllele = null;
	public Allele majorAllele = null;
	public Allele alternateAllele = null;
	public Allele refSeqAllele = null;
	public String ref,alt;
	
	//this is the percentage of samples in the samples list that both this SNP and its compared SNP had data for
	public float percentSamplesSharedWithCompSNP;
	
	public ArrayList<Allele> alleleList = new ArrayList<Allele>();	
	public Contig contig;	
	public int positionInCodon = -1;
	public float qualityScore = -1;
	
	//all of these are to enable filtering of SNPs with poor quality alt alleles (this happens when we have sequence specific errors)
	public int refAlleleQualSum;
	public int altAlleleQualSum;	
	public int refAlleleCount;
	public int altAlleleCount;
	public float meanQualAltAllele;
	public float meanQualRefAllele;
	//this is computed as Math.abs(((snp.meanQualRefAllele - snp.meanQualAltAllele)/snp.meanQualRefAllele)*100);
	//i.e. the percent difference between the average base qualities of the reference and alternate alleles
	public float percentDiffRefAlleleQualAltAlleleQual;
	
	//SNP effect, as Integer, see static constants in SNPEffectPredictor class
	public int snpEffect;
	//SNP effect, as String, provided by snpEff annotation
	public String snpEffAnotationString = null;	
	
	//the amino acid translation of this SNP in the context of a given frame, assuming the reference allele
	public String refAlleleTranslation;
	//the amino acid translation of this SNP in the context of a given frame, assuming the alternate allele
	public String altAlleleTranslation;
	
	public boolean isMultiAllelic;//has this SNP got more than 2 bona fide alleles?
	public boolean isMNP;//has this SNP got ref and alt alleles of equal length greater than 1?
	public boolean isIndel;//this is a hack - technically we should have a different class for this but this would require monumental code changes
	
	//flags for whether various filters have been passed
	public boolean passMAC ;//minor allele count
	public boolean passMAF ;//minor allele frequency
	public boolean passPctDiff ;//percent difference between base qualities of reference allele and alternate allele
	public boolean passQualScoreFilter ;//quality score filter
	public boolean passPctHetSamplesFilter; //% heterozygous samples in SNP filter - used for filtering out SNPs in homozygous organisms where heterozygosity is an indication of cross-mapping
	public boolean passPctFiltHetSamplesFilter; //% filtered heterozygous samples in SNP filter - used for filtering out SNPs in homozygous organisms where heterozygosity is an indication of cross-mapping
	public boolean passMultiAllelicFilter; //has this SNP got more than 2 bona fide alleles?
	public boolean passMonoMorphicAlt; //are all the samples homozygous for the alternate allele?
	
	//this evaluates to true when this SNP has at least one sample that also satisfies the global minor allele count AND frequency
	public boolean withinSampleCountsAndFreqFilterPassed = false;
	
	//counts for samples with major/minor allele as genotype
	public int numberOfSamplesWithMajorAlleleAsGenotype, numberOfSamplesWithMinorAlleleAsGenotype;
//	the same as percentages
	public float pctSamplesWithMajorAlleleAsGenotype, pctSamplesWithMinorAlleleAsGenotype;
	
	//a list of filter objects
	public ArrayList<SNPFilter> filters;
	
	//the percentage of samples for which this SNP contains no data
	public float pctMissingSamples;
	
	//the percentage of samples that have at least the read depth required by our filtering parameter 
	public float pctSamplesWithSufficientDepth;
	
	
	
	
	//===========================================================================================================================
	
	public void workOutMinorMajorAlleles()
	{
		
		
		//first find out how many alleles we have
		int numAlleles = alleleList.size();		
		//sort the list in descending order
		Collections.sort(alleleList,Collections.reverseOrder());
		
//		System.out.println("workOutMinorMajorAlleles");
//		for (Allele allele : alleleList)
//		{
//			System.out.println(allele + ": " + allele.count);
//		}
		
		//we can have one, two or more different alleles
		if(numAlleles == 1)
		{
			majorAllele = alleleList.get(0);
			minorAllele = null;
		}
		else if(numAlleles >= 2)
		{
			majorAllele = alleleList.get(0);
			minorAllele = alleleList.get(1);
		}
		
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public boolean isSameAsSNP(SNP compareSNP)
	{
		boolean match = this.refSeqName.equals(compareSNP.refSeqName) &&
						this.position == compareSNP.position;
		
		return match;
	}
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public String getPolymorphism()
	{
		return refSeqAllele + "/" + alternateAllele;
	}
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	@Override
	public String toString()
	{
		return "SNP at " + position + " on ref seq " + refSeqName;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//counts the alleles for String in this format:
	//G|2|80|+|2|80|-|0|0
	private static int countAlleles(String alleleCountsStr)
	{
		int indexFirstPipeChar = alleleCountsStr.indexOf("|");
		int indexSecondPipeChar = alleleCountsStr.indexOf("|",indexFirstPipeChar+1);
		String countStr = alleleCountsStr.substring(indexFirstPipeChar+1, indexSecondPipeChar);
		
		return Integer.parseInt(countStr);
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public boolean areHeterozygousSamplesPresent()
	{
		if(numHets == 0)
			return false;
		else
			return true;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public boolean isMonomorphicAlt()
	{

		boolean hom = true;
		
		//for this to be true, all sample genotypes need to be homozygous alternate allele
		//to check this we can simply look for zeros in their numeric genotype, and even if there is a single zero in one of these we 
		//can then return false
		for(SNPSample snpSample : samplesList)
		{
//			System.out.println("snpSample = " + snpSample);
//			System.out.println("snpSample.numericalGenotypeFromVCFString = " + snpSample.numericalGenotypeFromVCFString);
			
			if(snpSample.numericalGenotypeFromVCFString.contains("0"))
			{
				hom = false;
				break;
			}
		}
		
		
		return hom;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public float getMajorAlleleFrequency()
	{
		//divide the minor count by the total to get the ratio
		float freq = 0;
		if(majorAllele != null)
			freq = majorAllele.count / (float)totalReadCount;
				
		return freq;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public float getMinorAlleleFrequency()
	{
		//divide the minor count by the total to get the ratio
		float freq = 0;
		if(minorAllele != null)
			freq = minorAllele.count / (float)totalReadCount;

		return freq;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public SNPSample getSNPSampleByName(String name)
	{
		SNPSample sample = null;
		
		for (SNPSample snpSample : samplesList)
		{
			if(snpSample.name.trim().equals(name.trim()))
			{
				sample = snpSample;
				break;
			}
		}
		
		return sample;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
	public int computeNumberOfSamplesWithAlleleAsGenotype(Allele allele, int ploidy)
	{
		//the number of homozygous samples in this SNP that have this  allele  as their genotype
		int numSamplesWithAlleleAsGenotype = 0;		
		for(SNPSample sample : samplesList)
		{
			String refOrAltGenotypeInFullNotation = SNPUtils.getRefOrAltGenotypeInFullNotation(allele, ploidy);
			
//			System.out.println("\nsample = " + sample.name);
//			System.out.println("sample.genotypeFromVCFAsString = " + sample.genotypeFromVCFAsString);
//			System.out.println("refOrAltGenotypeInFullNotation = " + refOrAltGenotypeInFullNotation);
			
			if(sample.vcfGenotypeIsHomozygous() && sample.genotypeFromVCFAsString.equals(refOrAltGenotypeInFullNotation))
			{
//				System.out.println("incrementing count");
				numSamplesWithAlleleAsGenotype++;
			}		
		}

		return numSamplesWithAlleleAsGenotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
}
