# Commands
This scipt details the commands used to find genomic intervals with continuous read coverage and then define transcript boundaries 

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
- [Run operon prediction pipeline in R](#run-operon-prediction-pipeline-in-r)
	- [Create table containing coverage intervals](#create-table-containing-coverage-intervals)
	- [Create table containing read coordinates](#create-table-containing-read-coordinates)
	- [Add a column to reads that tabulates start and stop frequencies](#add-a-column-to-reads-that-tabulates-start-and-stop-frequencies)
	- [Assign reads to each coverage interval](#assign-reads-to-each-coverage-interval)
	- [Specify region for consideration and separate all read rows into data frames by coverage interval coordinate](#specify-region-for-consideration-and-separate-all-read-rows-into-data-frames-by-coverage-interval-coordinate)
	- [Simulate randomly distributed reads from each coverage region](#simulate-randomly-distributed-reads-from-each-coverage-region)
	- [Determine the number of high-frequency starts and ends for each coverage interval](#determine-the-number-of-high-frequency-starts-and-ends-for-each-coverage-interval)

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
# Run operon prediction pipeline in R
## Create table containing coverage intervals
```{R, eval = F}
k12lbgenefint=read.table("C:/Users/jarrett.lebov/OneDrive - University of Maryland School of Medicine/Documents/DirectRNAseq/OperonPrediction/f_genes_covinter.txt", header=T, sep = "\t")
k12lbgenerint=read.table("C:/Users/jarrett.lebov/OneDrive - University of Maryland School of Medicine/Documents/DirectRNAseq/OperonPrediction/r_genes_covinter.txt", header=T, sep = "\t")
finter=k12lbgenefint
rinter=k12lbgenerint
inter=rbind(finter,rinter);head(inter)
inter$strand=c(rep("+",length(finter$starts)),rep("-",length(rinter$starts)));head(inter)
colnames(inter)=c("ileft","iright","istrand");head(inter)
inter$interlen=inter$iright-inter$ileft;head(inter)
inter$istart=ifelse(inter$istrand=="+",inter$ileft,inter$iright);head(inter)
inter$iend=ifelse(inter$istrand=="+",inter$iright,inter$ileft);head(inter)
```
## Create table containing read coordinates
```{R, eval = F}
k12lbgenereads=read.table("C:/Users/jarrett.lebov/OneDrive - University of Maryland School of Medicine/Documents/DirectRNAseq/K12_LB_midlog/SequencingResults/MappingStats/gene_reads.bed")
reads=k12lbgenereads
reads=reads[,c(2,3,6)];head(reads)
colnames(reads)=c("left","right","strand");head(reads)
reads$readlen=reads$right-reads$left;head(reads)
reads$start=ifelse(reads$strand=="+",reads$left,reads$right);head(reads)
reads$end=ifelse(reads$strand=="+",reads$right,reads$left);head(reads)
head(reads[reads$strand=="-",])
```
## Add a column to reads that tabulates start and stop frequencies
```{R, eval = F}
# Create a lookup table that stores the start postiion frequencies
startlook=as.data.frame(table(reads$start));head(startlook)
# Create a list of start elements to be looked up
elements5=as.character(reads$start);head(elements5)
# Assign start posistion frequencies to a variable
pos5=startlook$Freq;head(pos5)
# Look through all elements and link each element with its corresponding frequency from the lookup table and assign frequencies to a vector variable
names(pos5)=startlook$Var1
allfreq5=unname(pos5[elements5])
# Create a new data column of each start positions frequency
reads$startfreq=allfreq5;head(reads)
# Repeat this to create an endfreq column of reads
endlook=as.data.frame(table(reads$end))
elements3=as.character(reads$end)
pos3=endlook$Freq
names(pos3)=endlook$Var1
allfreq3=unname(pos3[elements3])
reads$endfreq=allfreq3;head(reads)
```
## Assign reads to each coverage interval
```{R, eval = F}
install.packages("data.table")
library(data.table)
# create a column in inter that contains coverage interval coordinates
inter$icoords=paste(inter$ileft,inter$iright, sep = "-");head(inter)
#separate inter into plus/minus strand tables
iplus=inter[inter$istrand=="+",]
iminus=inter[inter$istrand=="-",]
# separate reads into plus/minus strand tables
plus=reads[reads$strand=="+",];head(plus);length(plus$left)
minus=reads[reads$strand=="-",];head(minus)
# Add columns to plus/minus read tables that contain the coverage interval coordinates that overlap each read
setDT(plus)[iplus, rintcoords := icoords, on = .(right >= ileft, left <= iright)]
setDT(minus)[iminus, rintcoords := icoords, on = .(right >= ileft, left <= iright)]
# Add columns to plus/minus read tables that contain the coverage interval lengths of intervals that overlap each read
plus[iplus, rintlen := interlen, on = .(right >= ileft, left <= iright)];head(plus)
minus[iminus, rintlen := interlen, on = .(right >= ileft, left <= iright)];head(minus)
# combine both plus and minus tables back into new reads table
reads=rbind(plus,minus);head(reads)
```
## Specify region for consideration and separate all read rows into data frames by coverage interval coordinate
```{R, eval = F}
# create a table containing a subset of reads
lbound=20000
rbound=30000
subreads=reads[reads$left>=lbound&reads$right<=rbound&reads$strand=="+"&is.na(reads$rintcoords)==F,];head(subreads)
# determine all coordinate sets in the subreads rintcoords column
rintcoordlev=unique(subreads$rintcoords);rintcoordlev
# initialize empty list of data.frames
intreads=vector(mode="list", length = length(rintcoordlev))
# iterate through my coordinate set list (intcoordlev) and for each element in intcoordlev, append all rows containing that element to a dataframe within intreads
for(i in seq_along(rintcoordlev)){
  intreads[[i]]=subreads[subreads$rintcoords==rintcoordlev[i]]
};intreads
# rename list elements acoording to their coordinates interval
names(intreads)=rintcoordlev;head(intreads)
```
## Simulate randomly distributed reads from each coverage region
```{R, eval = F}
reps=100 # designate the number or simulation replicates
freqstat=vector(mode="list", length = length(rintcoordlev)) # initialize empty vector for mean and sd stats for each coverage interval
for(k in seq_along(rintcoordlev)){
  begsmax=numeric(reps) # initialize empty vector to be populated by the max begin frequency observed for a given coverage interval
  finsmax=numeric(reps) # initialize empty vector to be populated by the max finish frequency observed for a given coverage interval
  for(i in seq_along(1:reps)){
    randlen=sample(intreads[[k]]$readlen, size = length(intreads[[k]]$left),replace = F) # sample read lengths for a given coverage interval
    begs=sapply(max(intreads[[k]]$rintlen)-randlen, function(j) sample(seq(1,j), 1)) # randomize the location of the begginings of reads for within a given coverage interval
    fins=begs+randlen # calculate where the random reads end
    randreads=data.frame(randlen,begs,fins) # establish a data frame that contains length beggining and end information for an interval's random reads
    beglook=as.data.frame(table(randreads$begs)) # create a table that stores begs in one column and their frequencies in another column (this is the lookup table)
    elements5=as.character(randreads$begs) # create a list of "begs" lables to look up
    pos5=beglook$Freq # create a list of frequencies corresponding to each begs key
    names(pos5)=beglook$Var1 # add each begs key element as a column header (name) for each corresponding frequency
    allfreq5=unname(pos5[elements5]) # create a list that is compiled by reporting the frequency associated with each begs key element header
    randreads$begfreq=allfreq5 # create a new column in randreads that reports each randreads' frequency
    finlook=as.data.frame(table(randreads$fins)) # repeat last 6 steps to determine the end frequencies
    elements3=as.character(randreads$fins)
    pos3=finlook$Freq
    names(pos3)=finlook$Var1
    allfreq3=unname(pos3[elements3])
    randreads$finfreq=allfreq3
    begsmax[i]=max(randreads$begfreq) # append the max begining frequency for achieved for a given repetition over a coverage interval
    finsmax[i]=max(randreads$finfreq) # append the max ending frequency for achieved for a given repetition over a coverage interval
  }
  boundfreqmean=c(mean(begsmax),mean(finsmax)) # determine the mean for the max frequencies for the begining/ending of reads for a given interval
  boundfreqsd=c(sd(begsmax),sd(finsmax)) # determine the sd for the max frequencies for the begining/ending of reads for a given interval
  boundfreqstats=data.frame(boundfreqmean,boundfreqsd);boundfreqstats # create a dataframe containing the useful statistics for determined for each coverage interval
  names(boundfreqstats)=c("mean","sd");boundfreqstats # rename the columns
  boundfreqstats$boundary=c("start","end");boundfreqstats # add a new column inidcating whether each row corresponds to a beginning or ending position
  freqstat[[k]]=boundfreqstats # append each data frame to freqstat
}
names(freqstat)=rintcoordlev;freqstat # name each dataframe in freqstat according to its coverage interval
```
## Determine the number of high-frequency starts and ends for each coverage interval
```{R, eval = F}
# calculate the number of high frequency starts per interval
scool=sapply(names(intreads), function(i) length(unique(intreads[[i]][intreads[[i]]$startfreq>=freqstat[[i]][1,1]+freqstat[[i]][1,2],]$left)));head(scool)
# calculate the number of high frequency ends per interval
ecool=sapply(names(intreads), function(i) length(unique(intreads[[1]][intreads[[1]]$endfreq>freqstat[[1]][2,1]+freqstat[[1]][2,2],]$right)));head(ecool)
```
