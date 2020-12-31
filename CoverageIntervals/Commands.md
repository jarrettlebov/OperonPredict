# Commands
This scipt details the commands used to find genomic intervals with consistent read coverage

<!-- MarkdownTOC autolink="true" -->

- [Directories](#directories)
- [Perform basecalling with guppy_basecaller](#perform-basecalling-with-guppy_basecaller)
- [Map reads to reference E. coli K12 genome](#map-reads-to-reference-e-coli-k12-genome)
- [Prepare reads for visualization on IGV](#prepare-reads-for-visualization-on-igv)
	- [Sort and index bam file](#sort-and-index-bam-file)
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
## Isolate forward and reverse strand reads
```{bash, eval = F}
cd "$Align_dir"
samtools view -F 4 -F 16 -o f_sorted_mapped_K12_directRNA.bam sorted_mapped_K12_directRNA.bam
samtools view -f 16 -o r_sorted_mapped_K12_directRNA.bam sorted_mapped_K12_directRNA.bam
```
## Index reads
```{bash, eval = F}
cd "$Align_dir"
samtools index f_sorted_mapped_K12_directRNA.bam
samtools index r_sorted_mapped_K12_directRNA.bam
```
# Determine coverage at each position in the genome
```{bash, eval = F}
cd "$Align_dir"
samtools depth -aa f_sorted_mapped_K12_directRNA.bam > f_coverage.txt
samtools depth -aa r_sorted_mapped_K12_directRNA.bam > r_coverage.txt
```
## Determine coverage at all non-zero positions (this might require a grep -v "0")
```{bash, eval = F}
samtools depth f_sorted_mapped_K12_directRNA.bam > f_non0coverage.txt
samtools depth r_sorted_mapped_K12_directRNA.bam > r_non0coverage.txt
```
## Isolate only the positions column
```{bash, eval = F}
cd $Align_dir
cut -f2 f_non0coverage.txt > f_non0pos.txt
cut -f2 r_non0coverage.txt > r_non0pos.txt
```
# Define intervals with coverage
## Invoke python script which isolates intervals with read coverage and excludes regions without coverage
```{bash, eval = F}
python $Scripts/non0CovInter.py f_non0pos.txt > f_covinter.txt
python $Scripts/non0CovInter.py r_non0pos.txt > r_covinter.txt
```

