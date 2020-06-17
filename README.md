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

## STEP 4 - Run NanoCicle
### Classify the soft-clipped read supporting chimeric eccDNA and soft-clipped supporting simple eccDNA
~~~bash
python NanoCircle_arg.py Classify -i BC09_hg19.bam
~~~
Which will be saved in a folder temp_reads containing both simple and complex reads in .bam format. 
### Create a .bai index for the read .bam
~~~bash
samtools index temp_reads/*.bam
~~~
### Identify Simple eccDNA using the coverage file and classified reads
~~~bash
python NanoCircle_arg.py Simple -i barcode_1000_cov.bed --ibam temp_reads/Simple_reads.bam -q 60 -o barcode_circles.bed
~~~
### Classify the soft-clipped read supporting chimeric eccDNA and soft-clipped supporting simple eccDNA
## STEP 5 - After Analysis
### Jaccard Index
calculating jaccard index for each individual circle compared to the estimated region with coverage
~~~bash
bedtools intersect -wao -a barcode_Simple_circles_1000.bed -b barcode_1000_cov.bed | head -10 | awk -v OFS='\t' '{print $1,$2,$3,($4/((($3-$2)+($10-$9))-$4))}'
~~~
