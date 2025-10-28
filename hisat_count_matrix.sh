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

mkdir -p hisat2 

# 1. Build HISAT2 index
hisat2-build '/path/to/GCF_000146045.2_R64_genomic_copy.fna' hisat2/genome_index

# 2. Align reads (example shown for all 4 replicates)
hisat2 -x /hisat2/genome_index -U conA_rep1_trimmed.fq -S Con_A1.sam
hisat2 -x /hisat2/genome_index -U conA_rep2_trimmed.fq -S Con_A2.sam
hisat2 -x /hisat2/genome_index -U conB_rep1_trimmed.fq -S Con_B1.sam
hisat2 -x /hisat2/genome_index -U conB_rep2_trimmed.fq -S Con_B2.sam

# 3. Convert SAM → BAM
samtools view -bS Con_A1.sam > Con_A1.bam
samtools view -bS Con_A2.sam > Con_A2.bam
samtools view -bS Con_B1.sam > Con_B1.bam
samtools view -bS Con_B2.sam > Con_B2.bam

# 4. Sort BAM files
samtools sort Con_A1.bam -o Con_A1_sorted.bam
samtools sort Con_A2.bam -o Con_A2_sorted.bam
samtools sort Con_B1.bam -o Con_B1_sorted.bam
samtools sort Con_B2.bam -o Con_B2_sorted.bam

# 5. Index BAM files
samtools index Con_A1_sorted.bam
samtools index Con_A2_sorted.bam
samtools index Con_B1_sorted.bam
samtools index Con_B2_sorted.bam

#6. Generating count matrix
featureCounts -a 'path/to/GCF_000146045.2_R64_genomic_copy.gtf' -o hisat_gene_counts.txt Con_A1_sorted.bam Con_A2_sorted.bam Con_B1_sorted.bam Con_B2_sorted.bam

