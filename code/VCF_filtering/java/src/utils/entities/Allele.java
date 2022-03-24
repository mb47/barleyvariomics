package utils.entities;

import java.util.*;

public class Allele implements Comparable
{
	//the value of this allele, e.g. "C"
	public char value;
	public String stringValue;
	
	public Allele(char value)
	{
		this.value = value;
		stringValue=new String(""+value);
	}
	
	public Allele(String alleleString)
	{
		if(alleleString.length() == 1)
			this.value = alleleString.toCharArray()[0];
		stringValue=alleleString;
	}
	

	//the count of the allele in a SNP
	public int count;
	
	@Override
	public String toString()
	{
		return stringValue;
	}

	@Override
	public int compareTo(Object o)
	{
		Allele allele2 = (Allele)o;
		
		if(count > allele2.count)
			return 1;
		else if(count < allele2.count)
			return -1;
		else	 //if(allele1.count == allele2.count)
			return 0;
	}
	
	public boolean equals(Allele anotherAllele)
	{
		if(value == anotherAllele.value)
			return true;
		else
			return false;
	}
}
