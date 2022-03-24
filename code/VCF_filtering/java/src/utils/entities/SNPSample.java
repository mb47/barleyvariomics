package utils.entities;

import java.util.*;
import utils.entities.*;


public class SNPSample
{
	public String name;
	
	public SNP snp;

	public Allele minorAllele = null;
	public Allele majorAllele = null;

	public Allele alternateAllele = null;
	public Allele refSeqAllele = null;
	
	public ArrayList<Allele> alleleList = new ArrayList<Allele>();
	
	public String genotypeFromVCFAsString;
	public String numericalGenotypeFromVCFString;
	public int singleDigitNumericalGenotypeFromVCFString;
	public float genotypeFromVCFPhredScore = -1;
	public String genotypeAlleleString = null;
	
	public int readCount;
	
	static final float withinSampleHeterozygosityThreshold = 0.10f;
//	static final int withinSampleHeterozygosityMinReadCount = 2;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	public int getMajorAlleleCount()
	{
		return majorAllele.count;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	public int getMinorAlleleCount()
	{
		if(minorAllele != null)
			return minorAllele.count;
		else 
			return 0;
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	public float getMinorAlleleFrequency()
	{
		return minorAllele.count / (float)(majorAllele.count + minorAllele.count);
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	public int getTotalReadCount()
	{
		return majorAllele.count + minorAllele.count;
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	public boolean isSameAsSample(SNPSample snpSample)
	{
		return name.equals(snpSample.name);
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	public String getGenotype()
	{
		String genotype = null;
		
		if(minorAllele != null)
			genotype = majorAllele + "/" + minorAllele;
		else if(minorAllele == null && majorAllele != null)
			genotype = majorAllele.toString();
		else
			genotype = "";
		
		//this is a hack to accommodate the samtools version of VCF
		if(genotypeAlleleString != null)
			genotype = genotypeAlleleString;
		
		return genotype;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//check for heterozygosity by counting the number of entities making up the VCF genotype call 
	public boolean vcfGenotypeIsHomozygous()
	{		
		//iterate over the components of the  raw VCF genotype String
		String[] numberTokens = numericalGenotypeFromVCFString.split("/");
		
		//this lookup holds the different entries - there should not be more than two (0 and/or 1, which is ref and alt respectively)
		//just put the 0 or 1 in as both key and value
		HashMap<String,String> lookup = new HashMap<String,String>();
		for (int i = 0; i < numberTokens.length; i++)
		{	
			String entry = numberTokens[i].trim();
			//check if this entry exists in the lookup			
			//if not, add it
			if(lookup.get(entry) == null)	
				lookup.put(entry, entry);
		}
		
		//check the size of the lookup -- if it only has one entry, then this sample is homozygous
		if(lookup.size() == 1)
			return true;
		//otherwise it's a het
		return false;
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public void workOutMinorMajorAlleles()
	{

		//first find out how many alleles we have
		int numAlleles = alleleList.size();		
		//sort the list in descending order
		Collections.sort(alleleList,Collections.reverseOrder());
			
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
	
//		System.out.println("minorAllele = " + minorAllele);
//		System.out.println("majorAllele = " + majorAllele);
		
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public boolean rawGenotypeIsHomozygous()
	{
		if(getMAFFilteredGenotype().contains("/"))
			return false;
		else
			return true;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
	public String getMAFFilteredGenotype()
	{
		String genotype = null;
		
		//if we have both alleles represented at a SNP in this sample then we need to call the genotype based on the expected error rate
		//to be conservative we should mandate that at least 10% of the reads within the sample contain the minor allele

		if(minorAllele != null)
		{
			//get the minor allele frequency for this SNP sample
			float minorAlleleFreqWithinSample = getMinorAlleleFrequency();
			//if it is below the threshold we call the genotype as major allele homozygous
			boolean minorAlleleFreqBelowCutoff = minorAlleleFreqWithinSample < withinSampleHeterozygosityThreshold;
			if(minorAlleleFreqBelowCutoff)
				genotype = majorAllele.toString();
			//otherwise we call a het
			else
				genotype = majorAllele + "/" + minorAllele;
		}
		//only major allele present
		else if(minorAllele == null && majorAllele != null)
			genotype = majorAllele.toString();
		//this will happen if this sample has no reads at this SNP, in which case we output nothing
		else
			genotype = "";

		
		return genotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
}
