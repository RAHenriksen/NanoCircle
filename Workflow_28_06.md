################################
###  ANALYSIS START 10-03-19 ###
###      FULL PIPELINE       ###
################################

########## STEP1
## adapter and barcode trimming 
~~~~bash
porechop -i bc05.reads.fastq -b bc05.barcode_trim -t 8
~~~~

#NB! husk lige fastq-stats hvis du vil have noget info omkring fastq filerne
assembly-stats -t */*.fastq > fastq_stats.txt

########## STEP2
## alignment
~~~~bash
minimap2 -t 8 -ax map-ont --secondary=no hg19.25chr.mmi BC05.fastq | samtools sort - > BC05.aln_hg19.bam
~~~~
minimap2 -t 6 -x map-ont -d GRCh38.mmi GRCh38.fa 
# -ax map-ont = Oxford Nanopore genomic reads
# --seconday=no With no reads mapped with SAM flag 0x100 (secondary flag). 
# hg19.25chr.mmi minimizer index for the reference


########## STEP3 - Trimming - before filtering
#convert fastq to fasta
sed '/^@/!d;s//>/;N' BC05.fastq > BC05.fasta
## python trimRollingCircle.py -t 5 -i $bc/chr1_14682287_14687166.fasta -o $bc/chr1_14682287_14687166.trim_rol.fasta
python /isdata/common/wql443/Thesis_2019/Analysis/trimRollingCircle.py -t 8 -i BC01.fasta -o BC01.trim.fasta

########## STEP4 - Identify eccDNA locus with mean depth >= 5
~~~~bash
bam=BC07.aln_hg19.bam
bname=$(basename $bam .aln_hg19.bam)
~~~~
### Convert bam to bed format and merge them as eccdna regions
~~~~bash
bedtools bamtobed -i $bam | bedtools sort | bedtools merge > $(basename $bam .bam.bed).merged.bed
~~~~
### Calculate mean of read depth in each eccdna region

bedtools genomecov -ibam BC07.aln_hg19.bam -d | awk 'BEGIN{OFS=FS="\t"}{if($3>0) print $1,$2-1,$2,$3}' > BC07.genomecov.bed

bedtools intersect -wo -a BC07.aln_hg19.merged.bed -b BC07.genomecov.bed | bedtools groupby -g 1,2,3 -c 7 -o mean > BC07.eccdna.bdg

### filtering out eccdna regions with mean depth >= 5 
~~~~bash
awk '$4>=5' $bname.eccdna.bdg > $(basename $bname.eccdna.bdg .bdg).ge_mean5.bdg
~~~~
awk '$4>=5' BC02.eccdna.bdg > BC02.ge_mean5.bdg
### extract reads from alignment result

~~~~bash
bc=BC01
for f in $bc; do awk -v f=$f '{print "samtools view "f".aln_hg19.bam "$1":"$2"-"$3" | cut -f1 | sort -u > "$1"_"$2"_"$3".txt ; python /isdata/common/wql443/Thesis_2019/Analysis/getFaLength.py -t filter -i "f".trim.fasta -incl "$1"_"$2"_"$3".txt -out fasta > "$1"_"$2"_"$3".fasta"}' $f.ge_mean5.bdg; done  > get_fasta_on_eccdna.sh


for f in $bc; do awk -v f=$f '{print "samtools view "f".aln_hg19.bam "$1":"$2"-"$3" | cut -f1 | sort -u > "$1"_"$2"_"$3".txt ; python /isdata/common/wql443/Thesis_2019/Analysis/getFaLength.py -t filter -i "f".trim.fasta -incl "$1"_"$2"_"$3".txt -out fasta > "$1"_"$2"_"$3".fasta"}' $f.ge_mean5.bdg; done  > get_fasta_on_eccdna.sh

bash get_fasta_on_eccdna.sh
~~~~

#then with these fasta files do the unicycler

########## STEP5 - ASSEMBLY
## assembly each eccdna region
# 1) unicycler
mkdir Unicycler 
for f in *.fasta; do unicycler -t 6 -l $f -o "$f"_uni; done

# mv succesful ones to this folder
find . -type f | cut -d/ -f2 | sort | uniq -d
mv $(find . -type f | cut -d/ -f2 | sort | uniq -d) Unicycler/
rm -rf *fasta_uni #as its just an empty log stating it failed
mv Unicycler/* .
for i in *.fasta *.txt; do mv "$i" "${i%.*}.fasta_uni/" ; done 
mv $(find . -type f | cut -d/ -f2 | sort | uniq -d) Unicycler/
# 2) Canu - Assembly
mkdir Canu
for f in *.fasta; do Gen_size=$(echo "$f" | awk -F'[_.]' '{print $3-$2}') ; canu -assemble -p assembly -d "$f"_canu_asm/ genomeSize=$Gen_size -nanopore-raw $f; done
# to see the number of files
find . -maxdepth 1 -type d -exec bash -c "echo -ne '{} '; ls '{}' | wc -l" \; |   awk '$NF==17' | wc -l

#to put fasta and txt file in to the canu folder - which allows for checking the size
for i in *.fasta *.txt; do mv "$i" "${i%.*}.fasta_canu_asm/" ; done 

#mv succes to Canu
mv $(find . -maxdepth 1 -type d -exec bash -c "echo -ne '{} '; ls '{}' | wc -l" \; |   awk '$NF==19') Canu/

# 3) Move the rest to failures, and check their circle size

mkdir Failures

########## CREATE A TABLE
# for all assembled + non assembled circles 
for f in *.fasta; do echo "$f" | awk -F'[_.]' '{print $1,$2,$3,$3-$2}'; done > Cirles_size.txt
ls *.fasta # in small screen to get one column with all coord 

### UNI
# see if they are contig or unitig
for f in *fasta_uni; do grep -e ">" $f/assembly.fasta ; done
# size of assembly
for f in *fasta_uni; do cat $f/assembly.gfa | awk -F':' '{print $3}' ; done

### CANU FULL
# see the length of the first assembly
for f in ./*canu_full; do head -1 $f/assembly.unitigs.fasta | awk -F'=' '{print $2}' ; done
# see if they are contig
for f in *canu_full; do grep -e ">" $f/assembly.unitigs.fasta | awk -F'=' '{print $2}' ; done
#for f in *canu_full; do cat $f/assembly.unitigs.fasta ; done

###### VALIDATION ######

#Comparison with some OF THE ASSEMBLYS BY COMPARING TO THE UNTRIMMED READS
# Print Read size
bioawk -c fastx '{print $name, length($seq)}' chr11_40022660_40027648.fasta | sort -k 2nr | awk '{print $2}'
# The 5 largest reads
bioawk -c fastx '{print $name, length($seq)}' chr11_40022660_40027648.fasta | sort -k 2nr | head -5
# Find read in a given region size, try similar size to the circle or the assembly
bioawk -c fastx '{print $name, length($seq)}' chr11_40022660_40027648.fasta | sort -k 2nr | awk '$2 < 5000 && $2 > 4000  {print ;}'
# Extract a read as fasta
samtools faidx chr11_40022660_40027648.fasta 98a6d695-a002-4bb8-84ed-274d99332add > chr11_40022660_40027648_read_98a6d.fa

samtools faidx BC05.fasta ffe33d83-0680-4ef2-bb46-eda948350f80 > chr11_40022660_40027648_read_98a6d.fa

# can do the same with trimmed reads

####### REMAPPING assembly to ref genome to se other potential mapping locations
# Denne er hurtigere 
minimap2 -t 6 -x map-ont --secondary=no hg19.25chr.mmi assembly.fasta

minimap2 -t 6 -x asm5 --secondary=no hg19.25chr.mmi assembly.fasta

# mmi file is an index file, and by remapping we're trying to explain how the assembled sequence make sense
minimap2 -t 6 -ax map-ont --secondary=no hg19.25chr.mmi assembly.fasta | samtools sort > assembly.bam
samtools index assembly.bam

### Remapping/Aligning of reads to the assembly??

### GENERAL HANDLING OF FASTA
#Print only read name
grep -o -E "^>\w+" chr2_82083496_82087081.trim_rol.fasta | tr -d ">"
# Create a text file with the good info
grep "^>" chr2_82083496_82087081.trim_rol.fasta | awk -v OFS=' ' '{print $1,$3,$7}' > chr2_82083496_82087081.read_info.txt

#Create a distribution
# Read size for a circle (could do it for all samples) also for trimmed reads to get a view of read distribution
bioawk -c fastx '{print $name, length($seq)}' chr11_40022660_40027648.fasta | sort -k 2nr | awk '{print $2}'


# Sorting the fasta or .fai 
sort -k2,2nr BC12.fasta.fai | head -5
sort -k2,2nr BC12.aln_hg19_trim.fasta | head -5

### Create ref genome for the circles
#creating a ref fasta in the region for dot plots. Copy a succes_coord file to the location and write the coordinates in the bed file
cat Success_coord.bed 
bedtools getfasta -fi ~/UAMS/hg19.fa -bed Success_coord.bed -fo Success_coord.fa
#then extract specific regions from the fasta
samtools faidx Success_coord.fa chr8:110063780-110068833 > chr8_110063780_110068833_ref.fa

samtools faidx Success_coord.fa chr1:243928620-243938331 > chr1:243928620-243938331_ref.fa


### FOR MY TEST SET ###
#Extract that given region
samtools view -b BC05.aln_hg19.bam "chr1:243928620-243938331" > chr1_243928620_243938331_region.bam

mkdir -p BC01/eccdna/assembly/{Unicycler,Canu,Failures}
for dir in */; do mkdir -- "*/{assembly,orig}"; done

###### ALIGNING USING MINIMAP2 #####
minimap2 -t2 --cs -x map-ont BC01.aln_hg19_test2.fasta chr11_40022660_40027648.fasta | paftools.js view - | less

minimap2 -t2 -x map-ont BC01.aln_hg19_test2.fasta chr11_40022660_40027648.fasta




######### Dette var det kode piroon gav for at lave bdg men af en eller anden grund kan den ikke klare man giver stdin som anden variabel i bedtools intersect
~~~~bash
bedtools genomecov -ibam $bam -d | awk 'BEGIN{OFS=FS="\t"}{if($3>0) print $1,$2-1,$2,$3}' | bedtools intersect -wo -a $bname.aln_hg19.merged.bed -b stdin | bedtools groupby -g 1,2,3 -c 7 -o mean > $bname.eccdna.bdg
~~~~
bedtools genomecov -ibam BC07.aln_hg19.bam -d | awk 'BEGIN{OFS=FS="\t"}{if($3>0) print $1,$2-1,$2,$3}' | bedtools intersect -wo -a BC07.aln_hg19.merged.bed -b stdin | bedtools groupby -g 1,2,3 -c 7 -o mean > BC07.eccdna.bdg

Samtools view -q 60 .bam | head -n 10 | gawk ‘{print $3,$4,$21}’ | less -S

#Her så tager ud baseret på flaget noget info omkring readet
samtools view -c -F0x4 XXX.sam 

#hvor F er antallet af mapped reads og ændre det til f så er det antallet af unmapped

#-f int Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with 0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with0' (i.e. /^0[0-7]+/) [0].

#-F INT Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with 0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with0' (i.e. /^0[0-7]+/) [0].

## BARE SE READ_ID OG CIGAR STRING FRA DEN ORIGINALE BAM FIL
samtools view chr1_243928620_243938331_region.bam | cut -f 1,6 > READ_CIGAR.txt
grep -F 'S' READ_CIGAR.txt > READ_SPLIT.txt

minimap2 -t 6 -ax map-ont --secondary=yes /isdata/common/wql443/Thesis_2019/GRCh38.fa BC09.fastq | samtools sort - > BC09_test.bam

# Find the mapped & unmapped sequence reads 
samtools view -b -f 4 BC09.hg19_m10.bam > unmapped.bam
samtools view -b -F 4 file.bam > mapped.bam

#strand 
samtools view -F 16 Reads.bam 
samtools view -f 16 Reads.bam

#
samtools view unmapped.bam | head -1 | awk '{print $1,$10}'

### BLASTALL
seqtk sample -s 10 BC09.fastq 10000 > BC09_rand_10000.fastq
#convert to fasta
blastall -p blastn -i BC09_rand_10000.fasta -d nt -o lol -a 10 -m 7
#then somehow use MEGAN

##### PYSAM LEG



#EXTRACT SINGLE READ FROM BAM
samtools view chr1_243928620_243938331_region.bam | grep 0d250a15-2e17-4fce-ad21-bda8e79d0457 | awk '{print length($10)}'

samtools view chr1_243928620_243938331_region.bam | grep 035d1b05-527e-4ff8-b4d9-cda3677347af | awk '{print length($10)}'

#Se de SA dataen
samtools view -h -f 2048 -o loooool.bam chr1_243928620_243938331_region.bam

samtools view -f 0x800


### MAP TEST
minimap2 -t 6 -x map-ont -d GRCh38.mmi /isdata/common/wql443/Thesis_2019/GRCh38.fa

minimap2 -t 6 -k 10 -x map-ont -d hg19_mini10.mmi hg19_ref2.fa

minimap2 -t 8 -k 10 -ax map-ont --secondary=yes hg19_mini10.mmi BC09.fastq | samtools sort - > BC09.hg19_m10_sec.bam

### DISTRIBUTION EXTRACTION
bioawk -c fastx '{print $name, length($seq)}' BC01.fasta | sort -k 2nr > BC01_reads.txt

bioawk -c fastx '{print $name, length($seq)}' BC01_trim.fasta | sort -k 2nr > BC01_trim_readsize.txt

cat BC01_trim_readsize.txt | awk '$2 < 1000' | wc -l

#number of mapped and unmapped reads
samtools view -c -F 4 BC01.aln_hg19.bam

samtools view -c -f 4 BC01.aln_hg19.bam






#PAF FORMAT
minimap2 -x ava-ont -t 8 hg19.fa /isdata/common/wql443/Thesis_2019/Analysis/2.8samples/BC01/BC01.fasta > /isdata/common/wql443/Thesis_2019/output.paf

#SIMULATION
read_analysis.py -i BC05.fasta -r hg19.25chr.mmi -a 'minimap2' -t 6

simulator.py circular -r sperm_plasmids/p4339.fasta -c training --max_len 50000 --min_len 150 -o Non_perfect
simulator.py circular -r sperm_plasmids/p4339.fasta -c training --max_len 50000 --min_len 150 --perfect
simulator.py circular -r sperm_plasmids/p4339.fasta -c training --max_len 50000 --min_len 150

#TRY TO EXTRACT REGION 
chr1:243928000-243939000
chr1:243928620-243938331
samtools faidx hg19.fa chr1:243928000-243939000 > chr1_243928000_243939000_circ.fa

samtools faidx hg19.fa chr1:chr1:243928620-243938331 > chr1_243928620_243938331_circ.fa

#jeg har jo allerede lavet et træning af mit reelle dataset
simulator.py circular -r chr1_243928000_243939000_circ.fa -c training -n 200 --max_len 15000 --min_len 200
#MAP det
minimap2 -t 8 -ax map-ont --secondary=no hg19.25chr.mmi Sim_read.fa | samtools sort - > Sim_circ.bam
#OVERVEJ LIGE AT INDSÆT COORD_SORT IGEN

#EXTRAHER FASTA READ LÆNGDE
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' Sim_read.fa

#kun for aligned -A1 henter også næste linje
 awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' Sim_read.fa | grep -A1 aligned | grep -v aligned


#FIND ANTALLET AF BASEPAIR DER VARIERER
bedtools intersect -wo -a Reg1000_1000.bed -b Simple_circles.bed -f 0.70 -r | awk '{print $1,$2,$3,$4,$5,$6,$7,$7-$6,$8-$4}' | sort -k1


#Håndtering af supp csv filer
sed -e 's/"//g' -e 's/;/\t/g' detected_eccDNAs_muscle.csv  > detected_eccDNAs_muscle.bed
awk -F'\t' '$14 == "hconf"' detected_eccDNAs_muscle.bed > muscle_hconf.bed

awk -v OFS='\t' '{print $2,$4,$5}' muscle_hconf.bed > muscle_hconf_3col.bed

bedtools merge -i muscle_lowqual_3col.bed -d 1000 > muscle_lowqual_merged.bed

sed -i '/^chrMT\b/d' muscle_hconf_merged.bed

#assembly w. illumina
#needs to extract the illumina data from the regions


for f in *.fasta; do unicycler -t 6 -1 /isdata/common/wql443/Thesis_2019/Illumina/Sample1_R1.fastq -2 /isdata/common/wql443/Thesis_2019/Illumina/Sample1_R2.fastq -l $f -o Unicycler/"$f"_uni; done

unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads.fastq.gz -o output_dir


~~~~bash

bwa mem -t 8 /isdata/common/wql443/Thesis_2019/ref_genome/bwa/hg19.fa /isdata/common/wql443/Thesis_2019/Illumina/Sample1_R1.fastq /isdata/common/wql443/Thesis_2019/Illumina/Sample1_R2.fastq > BC01.aln_hg19_ill.bam

samtools sort BC01.aln_hg19_ill.bam -o BC01.illumina.sort.bam

bc=BC01
for f in $bc; do awk -v f=$f '{print "samtools view "f".aln_hg19.bam "$1":"$2"-"$3" | cut -f1 | sort -u > "$1"_"$2"_"$3".txt ; python /isdata/common/wql443/Thesis_2019/Analysis/getFaLength.py -t filter -i "f".trim.fasta -incl "$1"_"$2"_"$3".txt -out fasta > "$1"_"$2"_"$3".fasta"}' $f.ge_mean5.bdg; done  > get_fasta_on_eccdna.sh



bc=BC01
for f in $bc; do awk -v f=$f '{print "samtools view "f".illumina.sort.bam "$1":"$2"-"$3" | cut -f1 | sort -u > "$1"_"$2"_"$3"_R1.txt ; python /isdata/common/wql443/Thesis_2019/Analysis/getFaLength.py -t filter -i /isdata/common/wql443/Thesis_2019/Illumina/trimmed/Sample1_R1.fastq -q -incl "$1"_"$2"_"$3"_R1.txt -out fastq > "$1"_"$2"_"$3"_R1.fastq"}' $f.ge_mean5.bdg; done  > get_fasta_on_eccdna_ill_R1.sh

bash get_fasta_on_eccdna_ill_R1.sh

for f in $bc; do awk -v f=$f '{print "samtools view "f".illumina.sort.bam "$1":"$2"-"$3" | cut -f1 | sort -u > "$1"_"$2"_"$3"_R2.txt ; python /isdata/common/wql443/Thesis_2019/Analysis/getFaLength.py -t filter -i /isdata/common/wql443/Thesis_2019/Illumina/trimmed/Sample1_R2.fastq -q -incl "$1"_"$2"_"$3"_R2.txt -out fastq > "$1"_"$2"_"$3"_R2.fastq"}' $f.ge_mean5.bdg; done  > get_fasta_on_eccdna_ill_R2.sh

#inde i den folder med on reads
for f in *.fasta; do unicycler -t 6 -1 ../Ill_reg_reads/${f%%.*}_R1.fastq -2 ../Ill_reg_reads/${f%%.*}_R2.fastq -l $f -o Unicycler/"$f"_uni; done




~~~~
${i%%.*}.mp3

bedtools intersect -f 0.9 -r -wa -wb -a Reg1000_20000.bed -b Simple_circles.bed | sort -k1,1 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$8-$4}' > Overlapping1000_20000.bed

bwa mem -t 8 /isdata/common/wql443/Thesis_2019/ref_genome/bwa/hg19.fa Sample1_R1.fastq Sample1_R2.fastq > BC01.aln_hg19_ill.bam


java -jar /isdata/common/wql443/Trimmomatic-0.39/trimmomatic-0.39.jar PE Sample1_R1.fastq Sample1_R2.fastq trimmed/Sample1_R1.fastq trimmed/Sample1_R1_unpaired.fastq trimmed/Sample1_R2.fastq trimmed/Sample1_R2_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36


cat */assembly.fasta | grep '>' |  wc -l
cat */assembly.fasta | grep '>' -A1 | wc -l

find */*assembly.fasta -type f -size 0


### For at lave en tilsvarende bdg fil for illumina dataen

bedtools genomecov -ibam BC02.illumina.sort.bam -d | awk 'BEGIN{OFS=FS="\t"}{if($3>0) print $1,$2-1,$2,$3}' > BC02.ill_cov.bed

bedtools bamtobed -i BC02.illumina.sort.bam | bedtools sort | bedtools merge > BC02.illumina.merged.bed


bedtools intersect -wo -a BC01.illumina.merged.bed -b BC01.ill_cov.bed | bedtools groupby -g 1,2,3 -c 7 -o mean | awk '$3-$2>=100' | awk '$4 >= 5' > BC01.illumina.100_mean5.bdg


#mapping Q
samtools view x.bam | awk '{$2 > 0}' | awk '{print $5}' 
#hvilket er forkert da flag = 0 bare siger den er unpaired men mapper til forward strand, so $2 != 4 (unmapped)
samtools view -c -F 4 BC01.aln_hg19.bam
#bare tag dem alle sammen



minimap2 -t 8 -ax map-ont --secondary=no /isdata/common/wql443/Thesis_2019/ref_genome/hg19.25chr.mmi assembly.fasta | samtools sort - | samtools view
