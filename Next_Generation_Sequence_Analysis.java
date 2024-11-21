package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

import mutationAnalysis.objects.MSA;
import mutationAnalysis.objects.EvolutionaryDistance;
import parser.CreateViralPileup;
import parser.ParseFastq;
import readAnalysis.Pileup2Consensus;
import readAnalysis.objects.Fastq;
import readAnalysis.objects.Sequence;
import readAnalysis.objects.TranslatedSequence;
import statistic.tools.ReadStats;
import translation.TranslateDNAsequence;
import externalTools.RunClustalw2;
import externalTools.RunPhylip;
import externalTools.RunPicard;
import externalTools.RunR;
import externalTools.RunSamtools;
import externalTools.objects.MpileupLineFiltered;
import fileHandling.StringBuilderWriter;
import externalTools.RunBWA;

/**
 * 
 ** DeepSeqAnalysis works on HDV ILLUMINA Paired Read NGS Data
 ** Filter the reads based on Quality, Removal of Adapter sequencing
 **	Creating sam -> bam -> consensus sequence 
 */

public class DeepSeqAnalysis {

	public static void main(String[] args) {

		Locale.setDefault(Locale.US);
		StringBuilder analysisInformation=new StringBuilder();
		String out="/home/users/zusman/HBV/Shotgun_sequencing_Patient5/";
		String args0= "/home/users/zusman/HBV/Shotgun_sequencing_Patient5/hbv_refseq.fasta";
		File reference=new File(args0);  //		File reference=new File(args[0]);
		File log=new File(out, "logs");
		if(!log.exists()){
			log.mkdir();
		}
		RunBWA.runBwaIndex(reference, out); // Creates index of reference file 
		String args1="/home/users/zusman/HBV/Shotgun_sequencing_Patient5/";
		File direc=new File(args1);
		File[] runs=direc.listFiles(); // list all fastq files 
		ArrayList<Sequence> all=new ArrayList<Sequence>();
		ArrayList<TranslatedSequence> allAA=new ArrayList<TranslatedSequence>();
		analysisInformation.append("Output folder for following analysis: "+out+"\n");
		analysisInformation.append("Reference sequence is "+reference.getAbsolutePath()+".\n");
		StringBuilder sb=new StringBuilder();
		StringBuilder outsb=new StringBuilder();
		int maxDepth=150000;
		int Phred=33;
		int Quality=20;
		int flagFilter=4;
		//double mmp=0.25;
		int qualThresholdPileupPosition=20;
		double minFreq=0.1;
		boolean ambi=false;
		boolean indel=false;
		int refLength=680;
		int minQual=20;
		//		ArrayList<Patient> exp=ParseExperimentFile.extractPatientInformation(new File(args[2]));


		for(File input1:runs){
			if(input1.getName().endsWith("R1.fastq")){ // make sure all paired reads files have R1.fastq and R2.fastq in end of name 
				//442-R1.fastq.
				String id1=input1.getName().split("[-\\.]")[0].trim();
				String id2=id1+"-R2";
				id1= id1+"-R1";

				String temp =input1.getPath();
				temp = temp.replace("R1.fastq", "R2.fastq");
				File input2 = new File(temp);

				System.out.println("Reading file 1");

				Fastq fq1=new Fastq(input1,id1);
				int readDepth=fq1.getReads().size();
				analysisInformation.append("Read depth for R1 read Prior to Filtering "+readDepth+".\n");
				System.out.println("Reading file 2");
				Fastq fq2=new Fastq(input2,id2);
				readDepth=fq2.getReads().size();
				analysisInformation.append("Read depth for R2 read Prior to Filtering "+readDepth+".\n");


				
				//PERFORM FILTERING Q=20

				analysisInformation.append("Minimum mean quality per read is set to "+Quality+".\n");

				analysisInformation.append("Trimming paired reads for "+ id1+","+id2+"\n");
				
				System.out.println("Trimming paired reads for "+ id1+","+id2);		
				
				
				//PERFORM TRIMMING USING TRIM GALORE

				Fastq.filterWithTrim_Galore(Phred,Quality,fq1,fq2,out);

				Fastq fq_filtered1=new Fastq(new File(out+"/"+id1+"_val_1.fq"), id1);
				int readDepth1=fq_filtered1.getReads().size();
				analysisInformation.append("Read Satistics for"+id1+"_R1 read after Filtering "+readDepth1+".\n");
				analysisInformation.append("Run\t NumberReads\t MeanLength\t MedianLength\t MinLength\t MaxLength\t SDLength\t MeanQuality\t MedianQuality\t MinQuality\t MaxQuality\t SDQuality\n");
				analysisInformation.append(fq_filtered1.statisticsToString());

				Fastq fq_filtered2=new Fastq(new File(out+"/"+id2+"_val_2.fq"), id2);
				int readDepth2=fq_filtered2.getReads().size();
				analysisInformation.append("Read Statistics for"+id2+"_R2 read after Filtering "+readDepth2+".\n");
				analysisInformation.append("Run\t NumberReads\t MeanLength\t MedianLength\t MinLength\t MaxLength\t SDLength\t MeanQuality\t MedianQuality\t MinQuality\t MaxQuality\t SDQuality\n");
				analysisInformation.append(fq_filtered2.statisticsToString());

				// CREATING ALIGNMENT
				//Creating a unique name for alignment file as both the pairs are used to make one alignment file
				// please note that file name must contain - before paired info i.e. in 442-R1_val.....

				String mode=RunBWA.ILLUMINA;
				String aliMode=RunBWA.UNIQUE;
				//	*/
				String UniqueReadName = "";
				if(fq_filtered1.getFastqFile().getName().contains("_val_1.fq")){
					String ReadId[] = fq_filtered1.getFastqFile().getName().split(Pattern.quote("-"));
					UniqueReadName = ReadId[0].trim(); //442
				}
				analysisInformation.append("Running BWA Aligner on reads "+UniqueReadName+" and reference "+reference.getName()+".\n");

			//	/*
				

				//CREATING SAM FILE
				File samFile=new File(out+"/"+UniqueReadName+"_align.sam"); // 442_align.sam
				System.out.println("Running Alignment for "+ UniqueReadName);
				RunBWA.runBwaAligner(reference,fq_filtered1, fq_filtered2, samFile, out);

				//CREATING BAM FILE 442_align.sam to  442_align.bam
				File bamFile=new File(out+"/"+UniqueReadName+"_align.bam"); // 442_align.bam
				System.out.println("Creating Bam File for "+ UniqueReadName);
				RunSamtools.sortNindex(bamFile, out);

				
				RunPicard.sortSam(samFile,bamFile,out); // sam is input file and bam is the converted output file -- out = logfolder
				analysisInformation.append("Running BWASort and BWAIndex on "+UniqueReadName+".\n");

				
				System.out.println("Using samfile to Markduplicate for "+ UniqueReadName);
				File MarkDuplicate_Outfile=new File(out+"/"+UniqueReadName+"_dedup.bam"); // 442_align.bam
				
				RunPicard.MarkDuplicates(bamFile,MarkDuplicate_Outfile,out); 
				System.out.println("Indexing Markduplicates for "+ UniqueReadName);
				
				RunPicard.BuildbamIndex(MarkDuplicate_Outfile,out);
				System.out.println("Process complete for "+ UniqueReadName);
				
				File MappedBam= new File (out+"/"+UniqueReadName+"_mapped.bam");
				RunSamtools.view_mapped(MarkDuplicate_Outfile, flagFilter, MappedBam, out);
				System.out.println("Indexing Mapped Bam File for "+ UniqueReadName);
				RunPicard.BuildbamIndex(MappedBam,out);
			//	/*

				//COVERAGE PLOT

				RunR.runCoveragePlot(RunSamtools.runMpileup(MappedBam, reference, out, maxDepth), new File("/home/users/zusman/rscripts/plotCoverageDistribution.r"), refLength, out);
		     	analysisInformation.append("Plotting coverage for "+UniqueReadName+" with R.\n");

	//*	  

		     	//CREATING CONCENSUS SEQUEUNCE
		     	
		     	ArrayList<MpileupLineFiltered> mpf=CreateViralPileup.viralPileup(RunSamtools.runMpileup(MappedBam, reference, out, maxDepth), qualThresholdPileupPosition);
				analysisInformation.append("Creating viral pileup for "+UniqueReadName+".\n");
				analysisInformation.append("Quality threshold is set to "+qualThresholdPileupPosition+".\n");

				sb=new StringBuilder();
				sb.append("sequence\t position\t reference\t mappedReads\t passedReads_"+qualThresholdPileupPosition+"\t alternatives\t absoluteQuantity\t relativeQuantity\n");
				for(MpileupLineFiltered m:mpf){
					sb.append(m.toString()+"\n");
				}
				StringBuilderWriter.write(new File(out, "pileup/"+UniqueReadName+".pileup"), sb);
				Sequence nucs=Pileup2Consensus.pileup2Consensus(mpf, minFreq, ambi, indel);
				analysisInformation.append("Construct consensus sequences from pileup for "+UniqueReadName+".\n");
				analysisInformation.append("Minimal frequency is set to "+minFreq+", ambiguity is set to "+ambi+" and indel calling is set to "+indel+".\n");
				sb=new StringBuilder();
				sb.append(nucs.toString());
				StringBuilderWriter.write(new File(out, "consensus/"+UniqueReadName+".cns"), sb);

				all.add(nucs);
				sb=new StringBuilder();

				int i=0;
				TranslatedSequence aa=null;
				while(true){
					if(i>2||(aa!=null&&!aa.getSequence().contains("*"))){
						break;
					}else{
						aa=TranslateDNAsequence.translate(nucs, i);
						++i;
					}
				}
				if(aa.getSequence().contains("*")){
					System.err.println("Could not find a reading frame for "+UniqueReadName+" without stops");
				}
				else{
					sb.append(aa.toString());
					StringBuilderWriter.write(new File(out,"consensus/"+UniqueReadName+"_aa.cns"), sb);
					analysisInformation.append("Translating "+UniqueReadName+" into protein sequence.\n");
					analysisInformation.append("Using reading frame "+aa.getReadingFrame()+".\n");
					allAA.add(aa);
				}

				outsb.append(fq_filtered1.statisticsToString());
	//	*/		//*///				
			}
		}

		StringBuilderWriter.write(new File(out+"logs/analysis.log"), analysisInformation);
	}

}