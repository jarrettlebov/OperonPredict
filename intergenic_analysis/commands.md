# Commands
Below is a set of commands used to explore reads mapping to intergenic regions of E. coli K12. 

<!-- MarkdownTOC autolink="true" -->

- [Directories](#directories)
- [Perform basecalling with guppy_basecaller](#perform-basecalling-with-guppy_basecaller)
- [Map reads to reference E. coli K12 genome](#map-reads-to-reference-e-coli-k12-genome)
- [Prepare reads for visualization on IGV](#prepare-reads-for-visualization-on-igv)
	- [Sort and index bam file](#sort-and-index-bam-file)
	- [Isolate forward and reverse strand reads](#isolate-forward-and-reverse-strand-reads)
	- [Index reads](#index-reads)
- [Create K12 non-t-rRNA ".bed" files with locus tag information](#create-k12-non-t-rrna-bed-files-with-locus-tag-information)
- [Create .bed files from *_gene_reads.bam files](#create-bed-files-from-_gene_readsbam-files)
- [Figure out which gene intervals over lap each read intervals and count the number of genes that over lap reads](#figure-out-which-gene-intervals-over-lap-each-read-intervals-and-count-the-number-of-genes-that-over-lap-reads)
- [Caclulate the number of genes that are overlapped by a read which also overlaps another gene](#caclulate-the-number-of-genes-that-are-overlapped-by-a-read-which-also-overlaps-another-gene)
	- [Number of genes each read over laps](#number-of-genes-each-read-over-laps)
	- [Number of forward reads overlapping more than 1 gene](#number-of-forward-reads-overlapping-more-than-1-gene)
	- [Forard strand](#forard-strand)
	- [Forward strand: number of genes overlaped by a read](#forward-strand-number-of-genes-overlaped-by-a-read)
	- [Reverse strand: number of genes overlaped by a read](#reverse-strand-number-of-genes-overlaped-by-a-read)
	- [Number of reverse reads overlapping more than 1 gene](#number-of-reverse-reads-overlapping-more-than-1-gene)
	- [Reverse strand](#reverse-strand)
- [Create files that contain the number of bp a read overlaps and the percentage of a gene a read covers](#create-files-that-contain-the-number-of-bp-a-read-overlaps-and-the-percentage-of-a-gene-a-read-covers)
	- [Produce a file with all forward bp overlap values in single column](#produce-a-file-with-all-forward-bp-overlap-values-in-single-column)
	- [Produce a file with all reverse bp overlap values in single column](#produce-a-file-with-all-reverse-bp-overlap-values-in-single-column)
	- [Produce a file with all forward and reverse bp overlap values in single column](#produce-a-file-with-all-forward-and-reverse-bp-overlap-values-in-single-column)
	- [Produce a file with all forward bp overlap values for multigene reads in single column](#produce-a-file-with-all-forward-bp-overlap-values-for-multigene-reads-in-single-column)
	- [Produce a file with all reverse bp overlap values for multigene reads in single column](#produce-a-file-with-all-reverse-bp-overlap-values-for-multigene-reads-in-single-column)
	- [Produce a file with all forward and reverse bp overlap values for multigene reads in single column](#produce-a-file-with-all-forward-and-reverse-bp-overlap-values-for-multigene-reads-in-single-column)
	- [Produce a file with all forward fraction overlap values in single column](#produce-a-file-with-all-forward-fraction-overlap-values-in-single-column)
	- [Produce a file with all reverse fraction overlap values in single column](#produce-a-file-with-all-reverse-fraction-overlap-values-in-single-column)
	- [Produce a file with all forward and reverse fraction overlap values in single column](#produce-a-file-with-all-forward-and-reverse-fraction-overlap-values-in-single-column)
	- [Produce a file with all forward fraction overlap values for multigene reads in single column](#produce-a-file-with-all-forward-fraction-overlap-values-for-multigene-reads-in-single-column)
	- [Produce a file with all reverse fraction overlap values for multigene reads in single column](#produce-a-file-with-all-reverse-fraction-overlap-values-for-multigene-reads-in-single-column)
	- [Produce a file with all forward and reverse fraction overlap values for multigene reads in single column](#produce-a-file-with-all-forward-and-reverse-fraction-overlap-values-for-multigene-reads-in-single-column)
- [Create a coverage plot files](#create-a-coverage-plot-files)
	- [Index reads](#index-reads-1)
	- [Create a depth file for all positions](#create-a-depth-file-for-all-positions)
	- [Average the depth of every 10 kb, 100 kb, 1000 kb](#average-the-depth-of-every-10-kb-100-kb-1000-kb)

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
cd $outdir
samtools sort mapped_K12_directRNA.bam -o sorted_mapped_K12_directRNA.bam
samtools index sorted_mapped_K12_directRNA.bam
```
## Isolate forward and reverse strand reads
```{bash, eval = F}
samtools view -F 4 -F 16 -o f_sorted_mapped_K12_directRNA.bam sorted_mapped_K12_directRNA.bam
samtools view -f 16 -o r_sorted_mapped_K12_directRNA.bam sorted_mapped_K12_directRNA.bam
```
## Index reads
```{bash, eval = F}
samtools index f_sorted_mapped_K12_directRNA.bam
samtools index r_sorted_mapped_K12_directRNA.bam
```
# Create K12 non-t-rRNA ".bed" files with locus tag information
```{bash, eval = F}
cd $Ref_dir
grep $'\t+' GCF_000005845.2_ASM584v2_genomic.gff | sed 's/;/\t/g' | grep -v "biotype=rRNA\|biotype=tRNA" | awk '{print $1"\t"$4"\t"$5"\t"$NF"\t"$7}' | grep "locus" | sort -k2 -n | uniq > $Align_dir/f_gene_feats.bed
grep $'\t-' GCF_000005845.2_ASM584v2_genomic.gff | sed 's/;/\t/g' | grep -v "biotype=rRNA\|biotype=tRNA" | awk '{print $1"\t"$4"\t"$5"\t"$NF"\t"$7}' | grep "locus" | sort -k2 -n | uniq > $Align_dir/r_gene_feats.bed
cd $Align_dir
awk '{print $1,"\t"$2,"\t"$3,"\t"$4}' f_gene_feats.bed > f_gene_feats.min.bed
awk '{print $1,"\t"$2,"\t"$3,"\t"$4}' r_gene_feats.bed > r_gene_feats.min.bed
```
# Create .bed files from *_gene_reads.bam files
```{bash, eval = F}
cd $Align_dir
samtools view -bL f_gene_feats.bed f_sorted_mapped_K12_directRNA.bam > f_gene_reads.bam
samtools view -bL r_gene_feats.bed r_sorted_mapped_K12_directRNA.bam > r_gene_reads.bam
bedtools bamtobed -i f_gene_reads.bam | sort -k2 -n | uniq > f_gene_reads.bed
bedtools bamtobed -i r_gene_reads.bam | sort -k2 -n | uniq > r_gene_reads.bed
awk '{print $1,"\t"$2,"\t"$3,"\t"$4}' f_gene_reads.bed > f_gene_reads.min.bed
awk '{print $1,"\t"$2,"\t"$3,"\t"$4}' r_gene_reads.bed > r_gene_reads.min.bed
```
# Figure out which gene intervals over lap each read intervals and count the number of genes that over lap reads
```{bash, eval = F}
cd $Align_dir
rm *gene_reads_overlap.txt
python $Scripts/bederlap_finder.py f_gene_feats.min.bed f_gene_reads.min.bed > f_gene_reads_overlap.txt
python $Scripts/bederlap_finder.py r_gene_feats.min.bed r_gene_reads.min.bed > r_gene_reads_overlap.txt
```
# Caclulate the number of genes that are overlapped by a read which also overlaps another gene
## Number of genes each read over laps
```{bash, eval = F}
sed 's/^.*features = //g' f_gene_reads_overlap.txt > f_num_gene-read_overlaps.txt
sed 's/^.*features = //g' r_gene_reads_overlap.txt > r_num_gene-read_overlaps.txt
cat f_num_gene-read_overlaps.txt r_num_gene-read_overlaps.txt > all_num_gene-read_overlaps.txt
```
## Number of forward reads overlapping more than 1 gene
```{bash, eval = F}
grep -v "features = 0\|features = 1" f_gene_reads_overlap.txt | wc
```
## Forard strand
```{bash, eval = F}
grep -v "features = 0\|features = 1" f_gene_reads_overlap.txt | sed "s/^.*, \['//g" | sed 's/].*//g' | tr -d ",'" | sed 's/ /\n/g' | sort -n | uniq | wc
```
## Forward strand: number of genes overlaped by a read
```{bash, eval = F}
sed 's/^.*locus_tag=//g' f_gene_reads_overlap.txt | awk '{print $1}' | sort -n | uniq | wc
```
## Reverse strand: number of genes overlaped by a read
```{bash, eval = F}
sed 's/^.*locus_tag=//g' r_gene_reads_overlap.txt | awk '{print $1}' | sort -n | uniq | wc
```
## Number of reverse reads overlapping more than 1 gene
```{bash, eval = F}
grep -v "features = 0\|features = 1" r_gene_reads_overlap.txt | wc
```
## Reverse strand
```{bash, eval = F}
grep -v "features = 0\|features = 1" r_gene_reads_overlap.txt | sed "s/^.*, \['//g" | sed 's/].*//g' | tr -d ",'" | sed 's/ /\n/g' | sort -n | uniq | wc
```
# Create files that contain the number of bp a read overlaps and the percentage of a gene a read covers
## Produce a file with all forward bp overlap values in single column
```{bash, eval = F}
sed 's/^.*locus_tag=//g' f_gene_reads_overlap.txt | sed 's/features =.*//g' | cut -d " " -f2- | sed 's/],.*//g' | tr -d "[," | sed 's/ /\n/g' > f_bp_gene-read_overlap.txt
```
## Produce a file with all reverse bp overlap values in single column
```{bash, eval = F}
sed 's/^.*locus_tag=//g' r_gene_reads_overlap.txt | sed 's/features =.*//g' | cut -d " " -f2- | sed 's/],.*//g' | tr -d "[," | sed 's/ /\n/g' > r_bp_gene-read_overlap.txt
```
## Produce a file with all forward and reverse bp overlap values in single column
```{bash, eval = F}
cat f_bp_gene-read_overlap.txt r_bp_gene-read_overlap.txt > all_bp_gene-read_overlap.txt
```
## Produce a file with all forward bp overlap values for multigene reads in single column
```{bash, eval = F}
grep -v "features = 0\|features = 1" f_gene_reads_overlap.txt | sed 's/^.*locus_tag=//g' | sed 's/features =.*//g' | cut -d " " -f2- | sed 's/],.*//g' | tr -d "[," | sed 's/ /\n/g' > f_bp_gene-read_multi-overlap.txt
```
## Produce a file with all reverse bp overlap values for multigene reads in single column
```{bash, eval = F}
grep -v "features = 0\|features = 1" r_gene_reads_overlap.txt | sed 's/^.*locus_tag=//g' | sed 's/features =.*//g' | cut -d " " -f2- | sed 's/],.*//g' | tr -d "[," | sed 's/ /\n/g' > r_bp_gene-read_multi-overlap.txt
```
## Produce a file with all forward and reverse bp overlap values for multigene reads in single column
```{bash, eval = F}
cat f_bp_gene-read_multi-overlap.txt r_bp_gene-read_multi-overlap.txt > all_bp_gene-read_multi-overlap.txt
```
## Produce a file with all forward fraction overlap values in single column
```{bash, eval = F}
sed 's/^.*locus_tag=//g' f_gene_reads_overlap.txt | sed 's/features =.*//g' | sed 's/^.*], //g' | tr -d "][," | sed 's/ /\n/g' | sed '/^$/d' > f_frac_gene-read_overlap.txt
```
## Produce a file with all reverse fraction overlap values in single column
```{bash, eval = F}
sed 's/^.*locus_tag=//g' r_gene_reads_overlap.txt | sed 's/features =.*//g' | sed 's/^.*], //g' | tr -d "][," | sed 's/ /\n/g' | sed '/^$/d' > r_frac_gene-read_overlap.txt
```
## Produce a file with all forward and reverse fraction overlap values in single column
```{bash, eval = F}
cat f_frac_gene-read_overlap.txt r_frac_gene-read_overlap.txt > all_frac_gene-read_overlap.txt
```
## Produce a file with all forward fraction overlap values for multigene reads in single column
```{bash, eval = F}
grep -v "features = 0\|features = 1" f_gene_reads_overlap.txt | sed 's/^.*locus_tag=//g' | sed 's/features =.*//g' | sed 's/^.*], //g' | tr -d "][," | sed 's/ /\n/g' | sed '/^$/d' > f_frac_gene-read_multi-overlap.txt
```
## Produce a file with all reverse fraction overlap values for multigene reads in single column
```{bash, eval = F}
grep -v "features = 0\|features = 1" r_gene_reads_overlap.txt | sed 's/^.*locus_tag=//g' | sed 's/features =.*//g' | sed 's/^.*], //g' | tr -d "][," | sed 's/ /\n/g' | sed '/^$/d' > r_frac_gene-read_multi-overlap.txt
```
## Produce a file with all forward and reverse fraction overlap values for multigene reads in single column
```{bash, eval = F}
cat f_frac_gene-read_multi-overlap.txt r_frac_gene-read_multi-overlap.txt > all_frac_gene-read_multi-overlap.txt
```
# Create a coverage plot files
## Index reads
```{bash, eval = F}
cd $Align_dir
samtools index f_gene_reads.bam
samtools index r_gene_reads.bam
```
## Create a depth file for all positions
```{bash, eval = F}
cd $Ref_dir
samtools depth -aa -r NC_000913.3:1-4641652 $Align_dir/f_gene_reads.bam $Align_dir/r_gene_reads.bam | awk '{print $0"\t"$3+$4}' > gene_feats_read_depth.txt
```
## Average the depth of every 10 kb, 100 kb, 1000 kb 
```{bash, eval = F}
awk '{sum+=$5}(NR%10000==0){mean=sum/10000;print mean; sum=0;next}' gene_feats_read_depth.txt > 10kb_gene_feats_read_depth.txt
awk '{sum+=$5}(NR%100000==0){mean=sum/100000;print mean; sum=0;next}' gene_feats_read_depth.txt > 100kb_gene_feats_read_depth.txt
awk '{sum+=$5}(NR%1000000==0){mean=sum/1000000;print mean; sum=0;next}' gene_feats_read_depth.txt > 1000kb_gene_feats_read_depth.txt
```
