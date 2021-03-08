# Commands
This scipt details the commands used to find genomic intervals with consistent read coverage

<!-- MarkdownTOC autolink="true" -->

- [Directories](#directories)
- [Perform basecalling with guppy_basecaller](#perform-basecalling-with-guppy_basecaller)
- [Map reads to reference E. coli K12 genome](#map-reads-to-reference-e-coli-k12-genome)
- [Prepare reads for visualization on IGV](#prepare-reads-for-visualization-on-igv)
	- [Sort and index bam file](#sort-and-index-bam-file)
	- [Extract read coordinates that map to non-t/rRNA regions](#extract-read-coordinates-that-map-to-non-trrna-regions)
	- [Isolate forward and reverse strand reads](#isolate-forward-and-reverse-strand-reads)
	- [Index reads](#index-reads)
- [Determine coverage at each position in the genome](#determine-coverage-at-each-position-in-the-genome)
	- [Determine coverage at all non-zero positions \(this might require a grep -v "0"\)](#determine-coverage-at-all-non-zero-positions-this-might-require-a-grep--v-0)
	- [Isolate only the positions column](#isolate-only-the-positions-column)
- [Define intervals with coverage](#define-intervals-with-coverage)
	- [Invoke python script which isolates intervals with read coverage and excludes regions without coverage](#invoke-python-script-which-isolates-intervals-with-read-coverage-and-excludes-regions-without-coverage)

<!-- /MarkdownTOC -->

# Directories
```{bash, eval = F}
Working_dir=/local/aberdeen2rw/julie/Jarrett/DirectRNAseq
Align_dir=/local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Align_out2
Ref_dir=/local/aberdeen2rw/julie/Jarrett/Ref_genomes
Programs=/local/aberdeen2rw/julie/Jarrett/Programs
Scripts=/local/aberdeen2rw/julie/Jarrett/Scripts
DB_dir=/local/aberdeen2rw/julie/Jarrett/Ref_genomes/blastdbs
FeatureCounts=/local/aberdeen2rw/julie/Jarrett/Programs/Subread/subread-2.0.1-source/bin/featureCounts
Miniconda3_bin_dir=/local/aberdeen2rw/julie/Jarrett/Programs/miniconda3/bin
Tombo_bin_dir=/usr/local/packages/python-3.5.2/bin/
Fast5_dir=/local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Jarrett/20200702_E_coli_K12_direct_totRNA/20200702_E_coli_K12_direct_totRNA/20200702_2003_MN21969_FAN42967_b833baa7/fast5/
fastqdir=/local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Jarrett/20200702_E_coli_K12_direct_totRNA/20200702_E_coli_K12_direct_totRNA/20200702_2003_MN21969_FAN42967_b833baa7/fastq2/pass
Test=/local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Jarrett/20200702_E_coli_K12_direct_totRNA/20200702_E_coli_K12_direct_totRNA/20200702_2003_MN21969_FAN42967_b833baa7/fastq2/Test
F5=/local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Jarrett/20200702_E_coli_K12_direct_totRNA/20200702_E_coli_K12_direct_totRNA/20200702_2003_MN21969_FAN42967_b833baa7/fastq2/Test/Single_Fast5_test/0
```
# Perform basecalling with guppy_basecaller
```{bash, eval = F}
guppy_basecaller --compress_fastq -i /local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Jarrett/20200702_E_coli_K12_direct_totRNA/20200702_E_coli_K12_direct_totRNA/20200702_2003_MN21969_FAN42967_b833baa7/fast5/ -s /local/aberdeen2rw/julie/Jarrett/DirectRNAseq/Jarrett/20200702_E_coli_K12_direct_totRNA/20200702_E_coli_K12_direct_totRNA/20200702_2003_MN21969_FAN42967_b833baa7/fastq2/ -c rna_r9.4.1_70bps_hac.cfg --fast5_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 --num_callers 8 --cpu_threads_per_caller 4
```
# Map reads to reference E. coli K12 genome
```{bash, eval = F}
combine_fastq(){
zcat "$fastqdir"/*.gz > "$Align_dir"/comb.fastq
}
combine_fastq

align_fastq(){
minimap2 -ax map-ont -t 4 "$Ref_dir/GCF_000005845.2_ASM584v2_genomic.fna" "$Align_dir"/comb.fastq | samtools view -bho "$Align_dir"/mapped_K12_directRNA.bam
}
align_fastq
```
# Prepare reads for visualization on IGV
## Sort and index bam file 
```{bash, eval = F}
cd "$Align_dir"
samtools sort mapped_K12_directRNA.bam -o sorted_mapped_K12_directRNA.bam
samtools index sorted_mapped_K12_directRNA.bam
```
## Extract read coordinates that map to non-t/rRNA regions
```{bash, eval = F}
cd "$Ref_dir"
grep -v "biotype=rRNA\|biotype=tRNA" GFF_files/GCF_000005845.2_ASM584v2_genomic.gff | grep "CDS" | awk -F "\t" '{print $1,$4,$5,$7}' OFS="\t" | sort -k2 -n | uniq > all_gene_cds.bed
cd $Align_dir
samtools view -bL "$Ref_dir"/all_gene_cds.bed sorted_mapped_K12_directRNA.bam > gene_reads.bam
bedtools bamtobed -i gene_reads.bam > gene_reads.bed
```
## Isolate forward and reverse strand reads
```{bash, eval = F}
cd "$Align_dir"
samtools view -F 4 -F 16 -o f_gene_reads.bam gene_reads.bam
samtools view -f 16 -o r_gene_reads.bam gene_reads.bam
```
## Index reads
```{bash, eval = F}
cd "$Align_dir"
samtools index f_gene_reads.bam
samtools index r_gene_reads.bam
```
# Determine coverage at each position in the genome
```{bash, eval = F}
cd "$Align_dir"
samtools depth -aa f_gene_reads.bam > f_genes_coverage.txt
samtools depth -aa r_gene_reads.bam > r_genes_coverage.txt
```
## Determine coverage at all non-zero positions (this might require a grep -v "0")
```{bash, eval = F}
cd "$Align_dir"
samtools depth f_gene_reads.bam > f_genes_non0coverage.txt
samtools depth r_gene_reads.bam > r_genes_non0coverage.txt
```
## Isolate only the positions column
```{bash, eval = F}
cd $Align_dir
cut -f2 f_genes_non0coverage.txt > f_genes_non0pos.txt
cut -f2 r_genes_non0coverage.txt > r_genes_non0pos.txt
```
# Define intervals with coverage
## Invoke [python](https://github.com/jarrettlebov/OperonPredict/blob/main/CoverageIntervals/PythonScripts/non0CovInter.py) script which isolates intervals with read coverage and excludes regions without coverage
```{bash, eval = F}
cd $Align_dir
python $Scripts/non0CovInter.py f_genes_non0pos.txt > f_genes_covinter.txt
python $Scripts/non0CovInter.py r_genes_non0pos.txt > r_genes_covinter.txt
```

