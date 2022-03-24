package utils.snps;

import java.util.*;
import utils.entities.*;


public class SNPUtils
{
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static String replaceBaseWithN(String sequence, SNP snp, int snpPosition)
	{
		String newString = null;
		
		try
		{
			String replacement = "N";

			if(sequence.length() == snpPosition)
				newString = sequence.substring(0, snpPosition) + replacement;
			else
				newString = sequence.substring(0, snpPosition) + replacement + sequence.substring(snpPosition +1);
			
			//		System.out.println("replacement = " + replacement);
			//		System.out.println("snpPosition = " + snpPosition);
			//		System.out.println("sequence = " + sequence);
			//		System.out.println("sequence.length() = " + sequence.length());
			//		System.out.println("snpPosition = " + snpPosition);
			//		System.out.println("newString = " + newString);
		}
		catch (Exception e)
		{
			System.out.println("sequence = " + sequence);
			System.out.println("sequence.length() = " + sequence.length());
			System.out.println("snpPosition = " + snpPosition);
			
			e.printStackTrace();
		}
		
		return newString;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static String replaceBaseWithIUPACCode(String sequence, SNP snp, int snpPosition)
	{
//		System.out.println("\nreplaceBaseWithIUPACCode");		
		
		String newString = null;
		
		try
		{

			//first extract the polymorphism of the SNP
			String polymorphism = snp.getPolymorphism();
			
//		System.out.println("polymorphism = " + polymorphism);
			
			String replacement = null;
			
			//iterate over the IUPAC options
			//		R 	A or G
			//		Y 	C or T
			//		S 	G or C
			//		W 	A or T
			//		K 	G or T
			//		M 	A or C
			
			if((polymorphism.charAt(0) == 'G' && polymorphism.charAt(2) == 'A') || (polymorphism.charAt(0) == 'A' && polymorphism.charAt(2) == 'G'))
				replacement = "R";
			else if((polymorphism.charAt(0) == 'C' && polymorphism.charAt(2) == 'T') || (polymorphism.charAt(0) == 'T' && polymorphism.charAt(2) == 'C'))
				replacement = "Y";
			else if((polymorphism.charAt(0) == 'G' && polymorphism.charAt(2) == 'C') || (polymorphism.charAt(0) == 'C' && polymorphism.charAt(2) == 'G'))
				replacement = "S";
			else if((polymorphism.charAt(0) == 'T' && polymorphism.charAt(2) == 'A') || (polymorphism.charAt(0) == 'A' && polymorphism.charAt(2) == 'T'))
				replacement = "W";
			else if((polymorphism.charAt(0) == 'G' && polymorphism.charAt(2) == 'T') || (polymorphism.charAt(0) == 'T' && polymorphism.charAt(2) == 'G'))
				replacement = "K";
			else if((polymorphism.charAt(0) == 'C' && polymorphism.charAt(2) == 'A') || (polymorphism.charAt(0) == 'A' && polymorphism.charAt(2) == 'C'))
				replacement = "M";
			else if(polymorphism.charAt(0) == polymorphism.charAt(2))
				replacement = ""+polymorphism.charAt(0);
			
			//				System.out.println("replacement = " + replacement);
			//		System.out.println("snpPosition = " + snpPosition);
			//		System.out.println("sequence = " + sequence);
			
			//replace the base in the string and return it
//		System.out.println("sequence = " + sequence);
//		System.out.println("sequence.length() = " + sequence.length());
//		System.out.println("snpPosition = " + snpPosition);
			
			
			if(sequence.length() == snpPosition)
				newString = sequence.substring(0, snpPosition) + replacement;
			else
				newString = sequence.substring(0, snpPosition) + replacement + sequence.substring(snpPosition +1);
			
			//		System.out.println("newString = " + newString);
		}
		catch (Exception e)
		{
			System.out.println("sequence = " + sequence);
			System.out.println("sequence.length() = " + sequence.length());
			System.out.println("snpPosition = " + snpPosition);
			
			e.printStackTrace();
		}

		return newString;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//replaces polymorphism notation (e.g. A/T) with IUPAC ambiguity code (e.g. W)
	public static String replacePolymorphismWithIUPACCode(String polymorphism)
	{
		String replacement = null;
		try
		{
			//iterate over the IUPAC options
			//		R 	A or G
			//		Y 	C or T
			//		S 	G or C
			//		W 	A or T
			//		K 	G or T
			//		M 	A or C
			
			if((polymorphism.charAt(0) == 'G' && polymorphism.charAt(2) == 'A') || (polymorphism.charAt(0) == 'A' && polymorphism.charAt(2) == 'G'))
				replacement = "R";
			else if((polymorphism.charAt(0) == 'C' && polymorphism.charAt(2) == 'T') || (polymorphism.charAt(0) == 'T' && polymorphism.charAt(2) == 'C'))
				replacement = "Y";
			else if((polymorphism.charAt(0) == 'G' && polymorphism.charAt(2) == 'C') || (polymorphism.charAt(0) == 'C' && polymorphism.charAt(2) == 'G'))
				replacement = "S";
			else if((polymorphism.charAt(0) == 'T' && polymorphism.charAt(2) == 'A') || (polymorphism.charAt(0) == 'A' && polymorphism.charAt(2) == 'T'))
				replacement = "W";
			else if((polymorphism.charAt(0) == 'G' && polymorphism.charAt(2) == 'T') || (polymorphism.charAt(0) == 'T' && polymorphism.charAt(2) == 'G'))
				replacement = "K";
			else if((polymorphism.charAt(0) == 'C' && polymorphism.charAt(2) == 'A') || (polymorphism.charAt(0) == 'A' && polymorphism.charAt(2) == 'C'))
				replacement = "M";
			else if(polymorphism.charAt(0) == polymorphism.charAt(2))
				replacement = ""+polymorphism.charAt(0);
			
		}
		catch (Exception e)
		{
		
			e.printStackTrace();
		}

		return replacement;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static float calculateAlleleCountRatio(HashMap<Character, Integer> alleleCountsLookup)
	{
		float ratio = - 1;
		int total = getTotalAlleleCount(alleleCountsLookup);
		int minorCount = getMinorAlleleCount(alleleCountsLookup);
		
		//divide the minor count by the total to get the ratio
		ratio = minorCount / (float)total;
		
		return ratio;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static int getTotalAlleleCount(HashMap<Character, Integer> alleleCountsLookup)
	{
		int total = 0;
		
		for(Character allele : alleleCountsLookup.keySet())
			total += alleleCountsLookup.get(allele);
		
		return total;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
	public static int getMajorAlleleCount(HashMap<Character, Integer> alleleCountsLookup)
	{
		int majorCount = Integer.MIN_VALUE;
		
		//first work out the total number of alleles
		for(Character allele : alleleCountsLookup.keySet())
		{
			//get the count for this allele
			int count = alleleCountsLookup.get(allele);
			//compare this count to the minor count var
			if(count > majorCount)
				majorCount = count;
		}
		
		return majorCount;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
	public static Character getMajorAllele(HashMap<Character, Integer> alleleCountsLookup)
	{
		Character majorAllele = null;
		
		//first get the count for the minor allele
		int majorAlleleCount = getMajorAlleleCount(alleleCountsLookup);
		
		//iterate over the lookup and find the allele that matches this count
		//if the counts are the same we use the first allele we encounter as then it doesn't matter which one gets picked
		for (Character allele : alleleCountsLookup.keySet())
		{
			if(alleleCountsLookup.get(allele) == majorAlleleCount)
				majorAllele = allele;
		}
		
		return majorAllele;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
	public static Character getMinorAllele(HashMap<Character, Integer> alleleCountsLookup)
	{
		Character minorAllele = null;
		
		//first get the count for the minor allele
		int minorAlleleCount = getMinorAlleleCount(alleleCountsLookup);
		
		//iterate over the lookup and find the allele that matches this count
		//if the counts are the same we use the first allele we encounter as then it doesn't matter which one gets picked
		for (Character allele : alleleCountsLookup.keySet())
		{
			if(alleleCountsLookup.get(allele) == minorAlleleCount)
				minorAllele = allele;
		}
		
		return minorAllele;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
	public static int getMinorAlleleCount(HashMap<Character, Integer> alleleCountsLookup)
	{
		//		int minorCount = Integer.MAX_VALUE;
		//
		//		//first work out the total number of alleles
		//		for(Character allele : alleleCountsLookup.keySet())
		//		{
		//			//get the count for this allele
		//			int count = alleleCountsLookup.get(allele);
		//			//compare this count to the minor count var
		//			if(count < minorCount)
		//				minorCount = count;
		//		}
		
		ArrayList<Integer> counts = new ArrayList<Integer>(alleleCountsLookup.values());
		Collections.sort(counts, Collections.reverseOrder());	
		int minorCount;
		if(counts.size() >= 2)
			minorCount = counts.get(1);
		else
			minorCount = counts.get(0);
		
		return minorCount;
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//works out the genotype for this result
	public static String extractGenotype(HashMap<Character, Integer> alleleCountsLookup)
	{
		String genotype = "";
		
		int total = 0;
		int majorCount = Integer.MIN_VALUE;
		int minorCount = Integer.MAX_VALUE;
		Character majorAllele = null;
		Character minorAllele = null;
		
		//first work out the total number of alleles
		for(Character allele : alleleCountsLookup.keySet())
		{
			//get the count for this allele
			int count = alleleCountsLookup.get(allele);
			//			System.out.println("allele " + allele + ", count " + count);
			//add it to the total
			total += count;
			//compare this count to the minor count var
			if(count > majorCount)
			{
				majorCount = count;
				majorAllele = allele;
			}
			if(count < minorCount)
			{
				minorCount = count;
				minorAllele = allele;
			}
		}
		
		//		System.out.println("major allele: " + majorAllele);
		//		System.out.println("minor allele: " + minorAllele);
		
		if(minorCount == 0)
			genotype = majorAllele.toString();
		else
			genotype = majorAllele.toString() + "/" + minorAllele.toString();
		
		//		System.out.println("genotype = " + genotype);
		
		return genotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static String formatAsGenotype(ArrayList<Character> alleles)
	{
		String genotype = null;
		
		if(alleles.size() == 1)
			genotype = alleles.get(0).toString();
		else if(alleles.size() ==2)
			genotype = alleles.get(0).toString() + "/" + alleles.get(1).toString();
		else
			genotype = "insufficient data";
		
		return genotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static String complementSequence(String sequence)
	{
		StringBuilder sb = new StringBuilder();
		char [] chars = sequence.toCharArray();
		
		for (int i = 0; i < chars.length; i++)
		{
			sb.append(getComplementaryBase(chars[i]));
		}
		
		return sb.toString();
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static char getComplementaryBase(char base)
	{
		char comp = base;
		
		switch(base)
		{
			case 'A' : comp = 'T'; break;
			case 'C' : comp = 'G';break;
			case 'G' : comp = 'C';break;
			case 'T' : comp = 'A';break;
			
			//also return the non-base chars
			case '/' : comp = '/';break;
			case '[' : comp = '[';break;
			case ']' : comp = ']';break;
		}
		
		return comp;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static  String reverseSequence(String sequence)
	{
		int length = sequence.length();
		StringBuffer reverse = new StringBuffer(length);
		for (int i = (length - 1); i >= 0; i--)
			reverse.append(sequence.charAt(i));
		return reverse.toString();
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static  String reverseComplementSequence(String sequence)
	{
		//complement first
		String comp = complementSequence(sequence);
		
		//then reverse
		return reverseSequence(comp);
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static boolean doGenotypesMatch(String genotype1, String genotype2)
	{
		//first check whether this sequence contains a "/" character -- if yes we have a heterozygote
		boolean gt1IsHet = genotype1.contains("/");
		boolean gt2IsHet = genotype2.contains("/");
		boolean bothHets = gt1IsHet && gt2IsHet;
		boolean bothHoms = !gt1IsHet && !gt2IsHet;
		
		//		System.out.println("gt1IsHet = " + gt1IsHet);
		//		System.out.println("gt2IsHet = " + gt2IsHet);
		//		System.out.println("bothHets = " + bothHets);
		//		System.out.println("bothHoms = " + bothHoms);
		
		//if these are not the same we have a disagreement and return false
		if(!(bothHets || bothHoms))
			return false;
		else if(bothHets)
		{
			//the String we need to split is something like "A/G"
			char gt1allele1 = genotype1.charAt(0);
			char gt1allele2 = genotype1.charAt(2);
			char gt2allele1 = genotype2.charAt(0);
			char gt2allele2 = genotype2.charAt(2);
			
			boolean match = ((gt1allele1 == gt2allele1) && (gt1allele2 == gt2allele2)) ||
			((gt1allele1 == gt2allele2) && (gt1allele2 == gt2allele1));
			
			return match;
		}
		else if(bothHoms)
		{
			return genotype1.equals(genotype2);
		}
		
		return false;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static Map sortByValue(Map map)
	{
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator()
		{
			public int compare(Object o1, Object o2)
			{
				return ((Comparable) ((Map.Entry) (o1)).getValue()).compareTo(((Map.Entry) (o2)).getValue());
			}
		});
		
		Map result = new LinkedHashMap();
		for (Iterator it = list.iterator(); it.hasNext();)
		{
			Map.Entry entry = (Map.Entry) it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static Allele getAlleleByValue(char value, ArrayList<Allele> alleleList)
	{
		Allele allele = null;
		
		for(Allele al : alleleList)
		{
			if(al.value == value)
			{
				allele = al;
				break;
			}
		}
		
		return allele;
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static boolean doPolymorphismsMatch(String polymorphism1, String polymorphism2, boolean isReverseComplement)
	{
		boolean match = false;
		
		//the assumption is that the format of the polymorphisms string is e.g. "A/G"
		String allele1_1 = polymorphism1.substring(0, 1);
		String allele1_2 = polymorphism1.substring(2);
		String allele2_1 = polymorphism2.substring(0, 1);
		String allele2_2 = polymorphism2.substring(2);
		
		//we need to allow for potential swaps of allele1 and allele2 - no guarantee they are in the same order
		//we also need to check for reverse complement
		
		
		return match; 
	}

	//------------------------------------------------------------------------------------------------------------------------	

	//work out the actual genotype based on the VCF string in FORMAT sample info field
	//i.e. we want to translate the numerical genotype into a nucleid acid one
	//e.g. 0/1 becomes A/C etc. depending on the ref and alt allele
	public static String translateVCFGenotype(String gtString, SNP snp)
	{	
		String genotype = gtString;
	
		//replace all zeroes with the reference allele
		genotype = genotype.replaceAll("0", snp.vcfEntry.ref);
		//replace all ones with the alternate allele
		genotype = genotype.replaceAll("1", snp.vcfEntry.alt);
		
		return genotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------
	
	public static int convertNumericalVCFGenotypeToSingleDigit(String numericalVCFGenotype)
	{
		int singleDigitGenotype = 0;

		//split the genotype string, which could look like this 0/1 or 0/1/1/1
		String [] tokens = numericalVCFGenotype.split("/");
		for (int i = 0; i < tokens.length; i++)
		{
			//add each digit to the sum
			singleDigitGenotype += Integer.parseInt(tokens[i]);
		}
		
		return singleDigitGenotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	

	//work out genotype based on VCF string in FORMAT sample info field
	//return a numerical value
	// 0 = homozygous for reference allele
	//  >= 1 means het
	//actual value can vary with ploidy
	//e.g. tetraploid homozygous for alt allele is 4 
	//add 1 for each instance of an alternate allele in the SNP caller's genotype
	public static int getNumericalGenotype(String gtString)
	{
		int genotype = 0;
		
		// sum the numbers in the genotype string as provided by the SNP caller
		String[] numberTokens = gtString.split("/");
		for (int i = 0; i < numberTokens.length; i++)
		{	
			genotype += Integer.parseInt(numberTokens[i].trim());
		}
		
		return genotype;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	

	//loads the sample information into a suitable lookup where we can retrieve it
	//uses the format strings as keys, sample info as values
	public static HashMap<String, String> getSampleInfoLookup(String formatString, String sampleInfo)
	{
		//format: GT:GQ:DP:RO:QR:AO:QA:GL
		HashMap<String, String> lookup = new HashMap<String, String>();
		String [] formatTokens = formatString.split(":");
		String [] infoTokens = sampleInfo.split(":");
		for (int i = 0; i < formatTokens.length; i++)
			lookup.put(formatTokens[i], infoTokens[i]);
		
		return lookup;		
	}
	
	
	
	//------------------------------------------------------------------------------------------------------------------------

	//parse info string and put into map
	public static HashMap<String, String> parseInfoString(VCFEntry vcfEntry)
	{		
		HashMap<String, String> infoMap = new HashMap<String, String>();
		//info string looks like this:
		//AB=0.1;ABP=16.9077;AC=57;AF=0.59375;AN=96;AO=395;CIGAR=1X;DP=561;DPRA=1.68592;EPP=860.742;EPPR=363.475;HWE=-94.3876;LEN=1;MEANALT=1;MQM=255;MQMR=255;NS=48;NUMALT=1;ODDS=0.980829;PAIRED=0;PAIREDR=0;RO=166;RPP=860.742;RPPR=363.475;RUN=1;SAP=860.742;SRP=363.475;TYPE=snp;XAI=0;XAM=0;XAS=0;XRI=0;XRM=0.000847139;XRS=0.000847139;BVAR
		String [] keyValuePairs = vcfEntry.info.split(";");

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
				if (tokens[1].endsWith(";"))
					value = tokens[1].substring(0, (tokens[1].length()));
				else
					value = tokens[1];
				infoMap.put(tokens[0], value);
			}
		}
		
		return infoMap;
	}
	
	//------------------------------------------------------------------------------------------------------------------------	
	
	//turns a single allele (e.g. "A") into a String genotype representation (e.g. A/A or A/A/A/A for tetraploids)
	public static String getRefOrAltGenotypeInFullNotation(Allele allele, int ploidy)
	{
		StringBuffer genotype = new StringBuffer();
		
		//add the allele char to this string as many times as our ploidy value, separated by "/" chars
		for(int i = 0; i < ploidy; i++)
		{
			genotype.append(allele.stringValue);
			//only add the "/" if this is not the last iteration
			if(i < (ploidy - 1))
				genotype.append("/");
		}
		
		return genotype.toString();
	}

	//------------------------------------------------------------------------------------------------------------------------	
	
	
	
}
