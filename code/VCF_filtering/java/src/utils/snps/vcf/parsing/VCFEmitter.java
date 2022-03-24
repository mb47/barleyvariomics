package utils.snps.vcf.parsing;

import java.io.*;
import java.util.*;
import utils.entities.*;

public class VCFEmitter
{
	//===================================================================
	
	public static void outputWholeVCFFile(String [] headers, Collection<VCFEntry> entries, File outFile)
	{
		try
		{
			BufferedWriter writer  = new BufferedWriter(new FileWriter(outFile));
			
			outputVCFHeader(writer, headers);
			
			for (VCFEntry vcfEntry : entries)
				outputSingleVCFEntry(writer, vcfEntry);
			
			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	//===================================================================
	
	public static void outputVCFHeader(BufferedWriter writer, String [] headers) throws IOException
	{
		for (int i = 0; i < headers.length; i++)
			writer.write(headers[i] + "\t");
		
		writer.newLine();
	}
	
	//===================================================================
	
	public static void outputWholeVCFFile(LinkedList<String> doubleHashHeaders, String [] headers, Collection<VCFEntry> entries, File outFile)
	{
		try
		{
			BufferedWriter writer  = new BufferedWriter(new FileWriter(outFile));
			
			outputVCFHeader(writer, doubleHashHeaders, headers);
			
			for (VCFEntry vcfEntry : entries)
				outputSingleVCFEntry(writer, vcfEntry);
			
			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	//===================================================================
	
	public static void outputVCFHeader(BufferedWriter writer, LinkedList<String> doubleHashHeaders, String [] headers) throws IOException
	{
		//output the double hash headers
		for(String doubleHashHeader : doubleHashHeaders)
		{
			writer.write(doubleHashHeader);
			writer.newLine();
		}			
		
		//output the single hash headers
		for (int i = 0; i < headers.length; i++)
			writer.write(headers[i] + "\t");
		
		writer.newLine();
	}
	
	//===================================================================
	
	public static void outputVCFHeader(BufferedWriter writer, ArrayList<String> headers) throws IOException
	{
		for (String header : headers)
			writer.write(header);
	}
	
	
	//===================================================================
	
	public static void outputSingleVCFEntry(BufferedWriter writer, VCFEntry vcfEntry) throws IOException
	{
		writer.write(vcfEntry.contig + "\t");
		writer.write(vcfEntry.position + "\t");
		writer.write(vcfEntry.ID + "\t");
		writer.write(vcfEntry.ref + "\t");
		writer.write(vcfEntry.alt + "\t");
		writer.write(vcfEntry.qual + "\t");
		writer.write(vcfEntry.filter + "\t");
		writer.write(vcfEntry.info + "\t");
		writer.write(vcfEntry.format + "\t");
		for(String sampleInfo : vcfEntry.sampleInfos)
			writer.write(sampleInfo  + "\t");
		writer.newLine();
	}
	
	//===================================================================
	
}
