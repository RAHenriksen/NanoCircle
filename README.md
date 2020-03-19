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
minimap2 -t 6 -x map-ont -d GRCh38.mmi GRCh38.fa 
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
bedtools genomecov -bg -ibam BC09_39w/BC09_39.aln_hg19.bam | bedtools merge -d 1000 -i stdin | sort -V -k1,1 -k2,2n > BC09_39w/BC09_39_1000_cov.bed
~~~ 

