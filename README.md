# PollenSeq
 pipeline for phsing by single sperm sequencing

Weiyi Zhang and Arslan Tariq

2024-2-27

# 1. make project directories and sub-directories; Timming < 1 min
mkdir -p PollenSeq/{00_Softwares,01_Reads.data,\
    02_Genome.data,03_Alignment-sperm,03_Alignment-parent,\
    04_SNP.sperm,04_SNP.parent,parent_tmpdir,sperm_tmpdir,05_Final.SNP}
mkdir PollenSeq/01_Reads.data/{Raw_reads,Clean_reads,Clean_reads_parent, \
    Clean_reads_spermCell}
export $POLLEN_DIR=`realpath PollenSeq`

⚠️ Please keep in mind, the $POLLEN_DIR is your work directory.

# 2. Download data; Timing < 60 min
Downloading the test data for this pipeline using the following command:

wget https://figshare.com/ndownloader/files/44669011 -O ->> test_data.tar.gz
tar zxf test_data.tar.gz
mv test_data/test.genome.fa $POLLEN_DIR/02_Genome.data/
mv test_data/* $POLLEN_DIR/01_Reads.data/Raw_reads
The test dataset utilized in this study represents a subset of authentic tea single sperm sequencing data. We extracted chromosome 1 and 2 segments spanning 1-100 Mb from the tea assembly DASZ to serve as test_genome.fa. Within this dataset, we specifically selected 37 single sperm samples for analysis.

# 3. Build Docker container; Timing < 20 min
You can find a Dockerfile in the PolleSeq git repo (https://github.com/zwycooky/PollenSeq). Use this file and build a docker container.

docker build -t pollenseq-container .
⚠️ You can also manually install software and adding the executable software or scripts to your system's $PATH.

# 4. Adapter removing and quality control from the sequencing reads; Timing < 2 min
Removing the adapter and low quality reads using the following command:

#clean reads
docker run -v ${POLLEN_DIR}/01_Reads.data/Raw_reads/:/input/ \
    -v ${POLLEN_DIR}/PollenSeq/01_Reads.data/Clean_reads/:/output/ pollenseq-container \
    Clean_rawReads.pl /input /output 5 6
#move parent clean reads into Clean_reads_parent/
mv $POLLEN_DIR/PollenSeq/01_Reads.data/Clean_reads/FD.*.fq.gz \
    POLLEN_DIR/PollenSeq/01_Reads.data/Clean_reads/Clean_reads_parent
#move sperm clean reads into Clean_reads_spermCell/
mv $POLLEN_DIR/PollenSeq/01_Reads.data/Clean_reads/QSC.*.fq.gz \
    POLLEN_DIR/PollenSeq/01_Reads.data/Clean_reads/Clean_reads_spermCell
The number ‘5’ in the above command indicate threads using in fastp program, while the number ‘6’ indicate how many fastp will be running at the same time. You can adjust threads for your available cores. The clean reads of parent and sperm cell for this pipeline will be stored at $POLLEN_DIR/01_Reads.data/Clean_reads_parent/ and $POLLEN_DIR/01_Reads.data/Clean_reads_spermCell/, respectively.

# 5. Building index for genome; Timing < 5 min
To map the reads, we need to make an index of the reference genome. And for the SNP calling, we need dictionary files. The commands to make the index and dictionary are the following:

docker run -v ${POLLEN_DIR}/02_Genome.data/:/data pollenseq-container bwa index /data/test.genome.fa
docker run -v ${POLLEN_DIR}/02_Genome.data/:/data pollenseq-container samtools faidx /data/test.genome.fa
docker run -v ${POLLEN_DIR}/02_Genome.data/:/data pollenseq-container samtools dict -o /data/test.genome.dict /data/test.genome.fa
# 6. Mapping and SNP calling for Parent; Timing < 2h
First, reads of parent will be mapped to genome by bwa, then picard and samtools will be used to sort bam file, mark-duplicates, and indexing respectively. Second, SNPs will be called and filtered by GATK in parallel mode. SNP filtration is based on following parameters: QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0. Finally, split VCF files will be merged into final VCF file. Running the following command to perform mapping and SNP calling for parent:

#mapping parent reads
docker run -v ${POLLEN_DIR}/01_Reads.data/Clean_reads_parent/:/input/ \
    -v ${POLLEN_DIR}/02_Genome.data/:/genome/ \
    -v ${POLLEN_DIR}/03_Alignment-parent/:/output/ pollenseq-container 02_mapit.pl \
    -1 /input/FD.R1.clean.fq.gz \
    -2 /input/FD.R2.clean.fq.gz \
    -f /genome/test.genome.fa \
    -b /genome/test.genome.fa \ 
    -p 20 \
    -o /output/FD
#running SNP calling
docker run -v ${POLLEN_DIR}/03_Alignment-parent:/input/ \
  	-v ${POLLEN_DIR}/parent_tmpdir:/tmpdir/ \
  	-v ${POLLEN_DIR}/04_SNP.parent:/output/ \
  	-v ${POLLEN_DIR}/02_Genome.data/:/genome/ \
  	pollenseq-container SNP_caller.pl /input/ /tmpdir/ /output/ /genome/test.genome.fa 22
#merge SNP calling results
docker run -v ${POLLEN_DIR}/parent_tmpdir:/tmpdir/ \
  	-v ${POLLEN_DIR}/04_SNP.parent:/output/ \
	pollenseq-container merge_split_vcf.pl /tmpdir/PASS.vcf.list /output/Parent_merged.vcf
We perform SNP calling from parent’s mapping files using GATK HaplotypeCaller. Unfortunately, HaplotypeCaller is a single-threaded program. To run it parallel, we will divide our bam file into small chunks. SNP calling for Each contig or chromosome will be done in a sliding window method (10 Mb per window). All of these procedures can be automated using a script called SNP_caller.pl, eliminating the necessity to recreate the process.

⚠️ -p 20 and 22 is indicating 20 and 22 threads, respectively, you can change it according to your available cores.

# 7. Mapping and SNP calling for Sperm-sequencing data; Timing < 3.5 h
This process is to step 12 and also divided into multiple processes including read mapping, duplicate removal from mapping files, SNP calling, and filtering. SNP filtration is based on following parameters: QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0.

#mapping sperm cell reads
docker run -v ${POLLEN_DIR}/01_Reads.data/Clean_reads_spermCell/:/input/ \
    -v ${POLLEN_DIR}/02_Genome.data/:/genome/ \
    -v ${POLLEN_DIR}/03_Alignment-sperm/:/output/ pollenseq-container batch_mapit.pl \
    -i /input/ \
    -o /output/ \
    -f /genome/test.genome.fa \
    -p 5 \
    -t 6
#running SNP calling of sperm cell
docker run -v ${POLLEN_DIR}/03_Alignment-sperm:/input/ \
  	-v ${POLLEN_DIR}/sperm_tmpdir:/tmpdir/ \
  	-v ${POLLEN_DIR}/04_SNP.sperm:/output/ \
  	-v ${POLLEN_DIR}/02_Genome.data/:/genome/ \
  	pollenseq-container SNP_caller.pl /input/ /tmpdir/ /output/ /genome/test.genome.fa 22
#merge SNP calling results
docker run -v ${POLLEN_DIR}/sperm_tmpdir:/tmpdir/ \
  	-v ${POLLEN_DIR}/04_SNP.sperm:/output/ \
	pollenseq-container merge_split_vcf.pl /tmpdir/PASS.vcf.list /output/sperm_merged.vcf
⚠️ -p 5 indicates 5 threads are used in bwa mapping, while -t 6 represents running 6 bwa mapping at same time. 22 is indicating 22 threads for SNP calling, you can change it according to your available cores.

# 8. Second Round of SNP Filtration; Timing < 5 min
The accuracy of SNP calling of sperm cells directly affects the phasing results. Therefore, a stringent filter of sperm cell SNP calling is required. Specifically, only SNPs that are heterozygous in the parent are kept for further analysis, for this purpose SNPs called from the parent are used. Use the 01filter_snp.pl script from PollenSeq for this purpose. Samples with a relatively low SNP calling rate are considered as failing with whole genome amplification and need to be filtered out for further analysis. Use the stat_snps_V2.pl script from PollenSeq for this purpose. Finally, using script 02filter_snp.pl to keep SNPs showing heterozygous less than 2 samples for phasing, and heterozygous SNPs are setting to missing genotype (./.) in VCF. These steps can be executed using the following command:

#Only SNPs are heterozygous in the parent are kept
docker run -v ${POLLEN_DIR}/04_SNP.parent:/input/ \
    -v ${POLLEN_DIR}/04_SNP.sperm:/input2/ \
    -v ${POLLEN_DIR}/05_Final.SNP:/output/ \
    pollenseq-container bash -c "01filter_snp.pl /input/Parent_merged.vcf \
    /input2/sperm_merged.vcf > /output/sperm_cell_merged.01.vcf"
#filter out sample with SNP calling rate < 0.15
docker run -v ${POLLEN_DIR}/05_Final.SNP:/output/ pollenseq-container stat_snps_V2.pl \
    /output/sperm_cell_merged.01.vcf /output/Row_snp_stat.txt /output/Sample_snp_stat.txt
awk '$2 > 0.85' ${POLLEN_DIR}/05_Final.SNP/Sample_snp_stat.txt | cut -f1 | \
    tail -n +2 > ${POLLEN_DIR}/05_Final.SNP/snp_calling_rate_under_0.15.txt
#final filter
docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "02filter_snp.pl /input/sperm_cell_merged.01.vcf \
    /input/snp_calling_rate_under_0.15.txt > /input/sperm_cell_merged.02.vcf"
⚠️ Here, we use missed SNP calling rate of 0.85 as a threshold to filter out sperm samples with failing of whole genome amplification. You should set this threshold based on your own data. For the test data, no samples are filtered out.

# 9. split VCF files by contig/seq id; Timing < 1 min
This step is crucial for optimizing computational speed during the phasing process.

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "split_vcf_by_contig.pl /input/sperm_cell_merged.02.vcf \
    /input/split_vcf"
Now you have a directory named split_vcf in 05_Final.SNP.

# 10. Hapi-phasing for each VCF; Timing < 3 h
Running the Hapi-phasing step:

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container single.contig.pl /input/split_vcf /input/single_scripts /input/hapOut 20 2
⚠️ It will save the output of each job in the hapOut folder. Some contigs may not have a hap.txt file which can occur for 2 reasons: 1- Only one SNP was detected on this contig. 2- All SNPs are filtered out by the function hapiFilterError. 20 indicating 20 threads for Hapi phasing, while 2 indicating 2 Hapi phasing program are allowed running at the same time. You can adjust these parameters according to your available cores.

# 11. Crossover detection for each contig; Timing < 5 min
Detecting crossover events in each chromosome/contig by following command:

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ pollenseq-container single.co.pl /input/hapOut 2
⚠️ Crossovers detected in this step may not be the real crossover, misassembly or structure variation between reference genome and individual used for phasing will also lead to switches between the two haplotypes. Using the value "2" indicates the allocation of 2 threads for crossover detection. It's important to emphasize that the number of threads should not exceed the number of contigs or chromosomes to ensure meaningful analysis.

# 12. Running bin detection; Timing < 5 min
Calculating genetic bins in each chromosome/contig by following command:

#get chromosome/contig/scafforld length
docker run -v ${POLLEN_DIR}/02_Genome.data/:/data/ \
    -v ${POLLEN_DIR}/05_Final.SNP:/output/ \
    pollenseq-container bash -c "chromosome_length.pl -i /data/test.genome.fa \
    > /output/reference.length.txt"
#detect genetic bin
docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "find_bins.pl /input/hapOut \
    /input/reference.length.txt /input/results"
#merge binmap for each contig/chromosome
docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "merge_binmapMarker.pl /input/results/ \
    binmap.txt > /input/merged_binmap.txt"
#merge binMarkers for each contig/chromosome
docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "merge_binmapMarker.pl /input/results/ \
    binMarkers.txt > /input/merged_binMarkers.txt"
⚠️ Binmap file contain the position of each bin in the reference genome and binMarkers use ‘A’ and ‘B’ to refer the genotypes of bins, and ‘U’ indicates missing genotypes.

# 13. Running re-phasing bins; Timing < 1 min
When Hapi performs phasing in a misassembled region or between different contig/scaffolds, Hapi will randomly re-assign genotype, which will hamper long-range or chromosomal phasing. Therefore RephasingBins.R is designed to solve this issue by following command:

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "RephasingBins.R /input/merged_binMarkers.txt /input/"
⚠️ By default, a bin with missing genotypes > 30% will be removed before the re-phasing step. Re-phased bin markers will be saved in ‘merged_binMarkers.re-phasing.co8.missing30.txt’ which can be used for genetic map construction and bins which is re-phased will be saved in ‘wrong_phasing_bins.txt’.

# 14. Running genetic map construction; Timing < 1 min
Run the following command to construct genetic map:

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ \
    pollenseq-container bash -c "Mstmap.R \
    /input/merged_binMarkers.re-phasing.co8.missing30.txt /input/Linkage_group 0.000001"
It will save all the linkage groups in Linkage_group folder. We highly recommend you to carefully check the genetic map, especially in the situation where the physical map cannot match genetic map, which may be caused by a misassembly or wrong genotyping of bins.

⚠️ If it gives you error like this: Error in plot.new() : figure margins too large Calls: plot … image -> image.default -> plot -> plot.default -> plot.new Execution halted. Just ignore it. Or if you want the graphs then you can resolve this error by changing the margin of R plotting.
0.000001 indicating pvalue for grouping.

# 15. Correcting wrong-phased haplotypes; Timing < 1 min
Run the following command to correct the wrong-phased haplotypes:

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ pollenseq-container \
    bash -c "correct_haps.pl /input/merged_binMarkers.re-phasing.co8.missing30.txt \
    /input/merged_binmap.txt /input/wrong_phasing_bins.txt \
    /input/hapOut /input/Linkage_group /input/correct.haps.txt"
The output file ‘correct_haps.txt’ shows the chromosomal-level phased SNPs which can be used for ASE analysis, haplotype-resolved genome assembly, and other applications.

# 16. Scaffold bins to pseudomolecules-based genetic map; Timing < 1 min
The following command is used to scaffold bins to pseudomolecules based on genetic map:

docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ pollenseq-container \
    bash -c "Re-anchored.pl /input/Linkage_group \
    /input/merged_binmap.txt > /input/Re_anchored_genome.txt"
# 17. Running crossover detection in pseudomolecules; Timing < 1 min
Two steps are used to detect crossovers in pseudomolecules. First, linkage groups are used to identify crossover between two genetic bins. Secondly, crossover events are filtered on the basis of physical maps. For this purpose, the following commands should be entered:

#detect co
docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ pollenseq-container \
    bash -c "Detect_co_from_genetic_map.R \
    /input/merged_binMarkers.re-phasing.co8.missing30.txt \
    /input/Linkage_group /input/co_final.txt"
#co filteration
docker run -v ${POLLEN_DIR}/05_Final.SNP:/input/ pollenseq-container \
    bash -c "co_filter.R /input/co_final.txt \
    /input/Re_anchored_genome.txt /input/co_final.filter.txt"


