### NanoCircle, Rasmus Amund Henriksen, wql443@alumni.ku.dk

The github reporsitory for the under development
tool, NanoCircle 2020. Useful for identifying the coordinates
of both simple and chimeric circular molecules, sequenced
using long-read sequencing.

# Some presteps to perform before running NanoCircle

## STEP1 - Trimming and prehandling

### adapter and barcode trimming 
~~~~bash
porechop -i bc05.reads.fastq -b bc05.barcode_trim -t 8
~~~~
Use fastq-stats to obtain information regarding the sequences

## STEP2 - Alignment of sequence reads

### creating index
~~~~bash
minimap2 -t 6 -x map-ont -d GRCh37.mmi GRCh37.fa
~~~~

### Alignment
~~~bash
minimap2 -t 8 -ax map-ont --secondary=no hg19.25chr.mmi read_file.fastq | samtools sort - > barcode.aln_hg19.bam
# -ax map-ont = Oxford Nanopore genomic reads
# --seconday=no With no reads mapped with SAM flag 0x100 (secondary flag). 
# hg19.25chr.mmi minimizer index for the reference
~~~ 

## STEP3 - Identifying representative regions

### bedtools genomecov + merge
~~~bash
bedtools genomecov -bg -ibam barcode.aln_hg19.bam | bedtools merge -d 1000 -i stdin | sort -V -k1,1 -k2,2n > barcode_1000_cov.bed
~~~ 

# Running NanoCircle to identify the eccDNA coordinates
## STEP 4 - Classify the soft-clipped read supporting Simple eccDNA and soft-clipped supporting Chimeric eccDNA
~~~bash
python NanoCircle_arg.py Classify -i barcode_hg19.bam -d temp_reads
~~~
Which will be saved in a folder temp_reads containing both simple and complex reads in .bam format. 
### Create a .bai index for the read .bam
~~~bash
samtools index temp_reads/Simple_reads.bam
samtools index temp_reads/Chimeric_reads.bam
~~~
## STEP 5 - Identify Simple eccDNA using the coverage file and classified reads
~~~bash
python NanoCircle_arg.py Simple -i barcode_1000_cov.bed -b temp_reads/Simple_reads.bam -q 60 -o barcode_Simple_circles.bed
~~~
## STEP 6 - Identify Chimeric eccDNA using the coverage file and classified reads
~~~bash
python NanoCircle_arg.py Chimeric -i barcode_1000_cov.bed -b temp_reads/Chimeric_reads.bam -q 60 -o barcode_Chimeric_circles.bed
~~~
The output being a bed file with possible configurations of several chimeric eccDNA, since the identification extract reads originating from specific regions.
## STEP 7 - Merge Chimeric eccDNA configurations using the coverage file and classified reads
~~~bash
python NanoCircle_arg.py Merge -i barcode_Chimeric_circles.bed -o barcode_Merged_chimeric.bed
~~~
The output being a bed file with possible configurations of several chimeric eccDNA, since the identification extract reads originating from specific regions.
# Ideas not yet incorporated
## STEP 8 - Jaccard Index
calculating jaccard index for each individual circle compared to the estimated region with coverage
~~~bash
bedtools intersect -wao -a barcode_Simple_circles_1000.bed -b barcode_1000_cov.bed | head -10 | awk -v OFS='\t' '{print $1,$2,$3,($4/((($3-$2)+($10-$9))-$4))}'
~~~
To check if there might be a small region in between the coordinates without any coverage ?
Or just use mean coverage

# Different unix command useful for data preparation, analysis and test
~~~bash
#Removing reads aligning to contamination sources, while still keeping the bam format.
samtools view -h BC10.aln_hg19.bam |grep -v '>N'| grep -v '>A' |samtools view -Sbo BC10.bam -
# No reads.
cat BC07.fastq | awk '{print $1}' | grep '@' | sort | uniq | wc –l
# No of mapped
samtools view -F 0x4 BC07/BC07.aln_hg19.bam | cut -f 1 | sort | uniq | wc -l

~~~
