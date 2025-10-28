# getting fastqc files
fastqc conA_rep1.fq conA_rep2.fq conB_rep1.fq conB_rep2.fq

# trimming with trimmomatics
mkdir -p trimmed_reads
for file in conA_rep1.fq conA_rep2.fq conB_rep1.fq conB_rep2.fq; do
    java -jar /path/to/java/trimmomatic-0.39.jar SE \
      -phred33 \
      $file \
      ${file%.fq}_trimmed.fq \
      HEADCROP:12
    echo "Completed trimming for: $file"
done

# getting phred scores for trimmed files
fastqc conA_rep1_trimmed.fq conA_rep2_trimmed.fq conB_rep1_trimmed.fq conB_rep2_trimmed.fq

# indexing with bowtie2
mkdir -p bowtie2_index
bowtie2-build 'path/to/GCF_000146045.2_R64_genomic_copy.fna' bowtie2_index/genome_index

# aligning with reference genome
bowtie2 -x /bowtie2_index/genome_index -U conA_rep1_trimmed.fq -S conA_rep1.sam
bowtie2 -x /bowtie2_index/genome_index -U conA_rep2_trimmed.fq -S conA_rep2.sam
bowtie2 -x /bowtie2_index/genome_index -U conB_rep1_trimmed.fq -S conB_rep1.sam
bowtie2 -x /bowtie2_index/genome_index -U conB_rep2_trimmed.fq -S conB_rep2.sam

# converting to .bam files
samtools view -bS conA_rep1.sam > conA_rep1.bam
samtools view -bS conA_rep2.sam > conA_rep2.bam
samtools view -bS conB_rep1.sam > conB_rep1.bam
samtools view -bS conB_rep2.sam > conB_rep2.bam

# sorting bam files
samtools sort conA_rep1.bam -o conA_rep1_sorted.bam
samtools sort conA_rep2.bam -o conA_rep2_sorted.bam
samtools sort conB_rep1.bam -o conB_rep1_sorted.bam
samtools sort conB_rep2.bam -o conB_rep2_sorted.bam

# indexing bam files
samtools index conA_rep1_sorted.bam
samtools index conA_rep2_sorted.bam
samtools index conB_rep1_sorted.bam
samtools index conB_rep2_sorted.bam

# basic alignment stats for alignment quality checking
samtools flagstat conA_rep1_sorted.bam
samtools flagstat conA_rep2_sorted.bam
samtools flagstat conB_rep1_sorted.bam
samtools flagstat conB_rep2_sorted.bam

# per-chromosome counts
samtools idxstats conA_rep1_sorted.bam
samtools idxstats conA_rep2_sorted.bam
samtools idxstats conB_rep1_sorted.bam
samtools idxstats conB_rep2_sorted.bam

# counting mapped reads for each gene using the annotation file and sorted .bam files
featureCounts -a 'path/to/GCF_000146045.2_R64_genomic_copy.gtf' -o gene_counts.txt conA_rep1_sorted.bam conA_rep2_sorted.bam conB_rep1_sorted.bam conB_rep2_sorted.bam
