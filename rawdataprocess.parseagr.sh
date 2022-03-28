#!/bin/bash

#rawdatapath="/mnt/c/Users/ZJ/project/msb2017/rawdata/"
#tmp_dir="/mnt/c/Users/ZJ/project/msb2017/tmp_qc"
#samplelist="/mnt/c/Users/ZJ/project/msb2017/test.filereport_read_run_PRJNA348749_tsv.txt"
#workdir="/mnt/c/Users/ZJ/project/msb2017/results"
index="/mnt/c/Users/ZJ/project/msb2017/example_db/fasta/0x58v50"
gtffile="/mnt/c/Users/ZJ/project/msb2017/example_db/gff/0x58v50.gff"
#genomelist="/mnt/c/Users/ZJ/project/msb2017/ncbi-genomes-2022-03-22/GCF_000019425.1_ASM1942v1_genomic.len"

helpFunction()
{
   echo ""
   echo "Usage: $0 -s sample -1 read1.fastq.gz -2 read2.fastq.gz -o output_path -f filter_out_sense_and_antisense_reads -t num_of_threads"
   echo -e "\t-s one sample accession, eg. "
   echo -e "\t-1 path to read 1 of fastq files"
   echo -e "\t-2 path to read 2 of fastq files"
   echo -e "\t-o the directory to output results and temporary files"
   echo -e "\t-f to filter out sense and antisense reads of alignments or not, use true for yes and vice versa"
   echo -e "\t-t number of threads to generate alignments and some of the following steps"
   exit 1 # Exit script after printing help
}

while getopts "s:1:2:o:f:t:" opt
do
   case "$opt" in
      s ) sample="$OPTARG" ;;
      1 ) read1="$OPTARG" ;;
      2 ) read2="$OPTARG" ;;
      o ) workdir="$OPTARG" ;;
      f ) filterreads="$OPTARG" ;;
      t ) thread="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


if [ -z "$sample" ] || [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$workdir" ] || [ -z "$filterreads" ] || [ -z "$thread" ] 
then
   echo "Some or all of the parameters are empty";
   helpFunction
else

	echo -e "$sample"  | while read sample;
do
		if [[ $(ls $read1 | wc -l) -gt 0 ]] &&  [[ $(ls $read2 | wc -l) -gt 0 ]];
		then
			echo -e "start to run quality control on sample $sample at $(date)"
			mkdir $workdir/$sample
			mkdir $workdir/$sample/tmp_dir
			fastqc -o $workdir/$sample -f $read1 $read2 -d $workdir/$sample/tmp_dir -t $thread -f fastq
			echo -e "start to align reads of $sample to genome $index at $(date)"
			bwa mem -t $thread $index $read1 $read2 > $workdir/$sample/$sample.sam
			echo -e "start to process alignments on sample $sample at $(date)"
			samtools view -b -S -@ $thread $workdir/$sample/$sample.sam >  $workdir/$sample/$sample.bam
			cat $workdir/$sample/$sample.bam  | samtools sort -@ $thread -o $workdir/$sample/${sample}_sorted.bam
			samtools stats -@ $thread  $workdir/$sample/${sample}_sorted.bam >  $workdir/$sample/${sample}_sorted.stats
			samtools index -@ $thread $workdir/$sample/${sample}_sorted.bam $workdir/$sample/${sample}_sorted.bai
			#echo -e "complete raw data processing on $sample at $(date)"

			echo -e "start to summarize transcriptome profiles on sample $sample at $(date)"
		#bedtools bamtobed -i $workdir/$sample/${sample}_sorted.bam > $workdir/$sample/${sample}_sorted.bed
			
			if [[ $filterreads == "false" ]];
			then
				samtools depth  $workdir/$sample/${sample}_sorted.bam >  $workdir/$sample/${sample}_sorted.depth 
				htseq-count --format=bam --order=name --stranded=reverse --mode=union --type=gene --idattr=Name $workdir/$sample/${sample}_sorted.bam $gtffile > $workdir/$sample/$sample.htseq_report.txt
			elif  [[ $filterreads == "true" ]];
			then
				echo -e "spliting sense and antisense reads for sample $sample at $(date)"
				mkdir  $workdir/$sample/sense
				mkdir  $workdir/$sample/antisense

				samtools view -f 83,163  -b -S -@ $thread $workdir/$sample/${sample}.bam > $workdir/$sample/sense/${sample}_sense.bam
				samtools view -f 99,147 -b -S -@ $thread $workdir/$sample/${sample}.bam > $workdir/$sample/antisense/${sample}_antisense.bam
				echo -e "sense\nantisense" | while read sen
			do
				echo -e "start to generate transcription profile on ${sen} at $(date)"
				cat $workdir/$sample/${sen}/${sample}_${sen}.bam  | samtools sort -@ $thread -o $workdir/$sample/${sen}/${sample}_${sen}_sorted.bam
				samtools stats -@ $thread  $workdir/$sample/${sen}/${sample}_${sen}_sorted.bam >  $workdir/$sample/${sen}/${sample}_${sen}_sorted.stats
                        samtools index -@ $thread $workdir/$sample/${sen}/${sample}_${sen}_sorted.bam $workdir/$sample/${sen}/${sample}_${sen}_sorted.bai
			samtools depth  $workdir/$sample/${sen}/${sample}_${sen}_sorted.bam >  $workdir/$sample/${sen}/${sample}_${sen}_sorted.depth

			htseq-count --format=bam --order=name --stranded=reverse --mode=union --type=gene --idattr=Name $workdir/$sample/${sen}/${sample}_${sen}_sorted.bam $gtffile > $workdir/$sample/${sen}/${sample}_${sen}.htseq_report.txt

		done
		

			else
				echo -e "invalid arguement for -f filter_out_sense_and_antisense_reads, exiting analysis at $(date)"
			fi
			
			echo -e "transcriptomic analyses on sample $sample completed at $(date)"

			rm $workdir/$sample/$sample.sam
		#bedtools genomecov -bga -i $workdir/$sample/${sample}_sorted.bed  -g $genomelist  > $workdir/$sample/${sample}_cov.txt
	else
		echo -e "cannot find one or both fastq files for $sample, please check their paths for the command lines"
	fi

done

fi

