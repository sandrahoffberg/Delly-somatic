[![Code Ocean Logo](images/CO_logo_135x72.png)](http://codeocean.com/product)

<hr>

# Somatic structural variant & CNV calling in Delly

Delly is an integrated structural variant (SV) prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read and long-read massively parallel sequencing data. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome.  

This capsule does Structural Variant calling based on paired-end and split-reads and Copy Number Variation calling based on read depth in (Illumina) short reads from case (e.g., cancer) and matched control samples.

For more information, see the [Delly Github](https://github.com/dellytools/delly) page. 

Delly includes many tunable parameters (e.g., for discovering variants with minimum quality, within specific window sizes, etc). This capsule leaves most parameters as default.

## Input 
- BAM file that has been trimmed, aligned to reference, sorted, duplicate marked, read groups adjusted if necessary, indexed.  This capsule downloads data from an S3 bucket if a URL is provided in the App Panel. If no URL is provided, it will search for bam files in the **/data** directory.  It can use sample data from both locations in the same run. 
- In the **/data** directory, a genome in fasta format (```.fa``` ending) with an index created in Samtools that ends in ```.fa.fai```   is required to identify split-reads
- In the **/data** directory, a samplesheet, a tab-delimited (```.tsv```) sample description file where the first column is the sample ID (as in the VCF/BCF file) and the second column is either "tumor" or "control".
- In the **/data** directory, a comparesheet, a comma-delimited (```.csv```) file that has two samples on each line, the first is the control, the second is the matched tumor.  Sample names may appear on more than 1 line in the case of multiple cases matched to the same control or multiple controls matched to the same case. 
- Optional: In the **/data** directory, a mappability map, [downloaded here](https://gear.embl.de/data/delly/).  This is a ```fasta``` format file with repetitive regions of the genome hardmasked, and both an ```.fai``` and ```.gzi``` index. 
- Optional: In the **/data** directory, a BED file of repetitive regions to exclude (e.g., telomeres and centromeres).  The one used here was downloaded from [SpeedSeq Github page](https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed) and linked to from smoove Github page.

## App Panel Parameters

### Main Parameters

Path to samplesheet
- A tab-delimited (```.tsv```) sample description file where the first column is the sample ID (as in the VCF/BCF file) and the second column is either "tumor" or "control". [Default: samplesheet.tsv in /data directory]

The URL of the S3 bucket with ```bams```
- If provided, will pull .bam files from this S3 location. [Default: s3://codeocean-public-data/example_datasets/chr11_tumor_bams/]

S3 Directory
- The location to save the ```bam``` files from S3 [Default: /scratch/data]

Path to comparesheet
- A comma-delimited (.csv) file that has two samples on each line, the first is the control, the second is the matched tumor. 
Sample names may appear on more than 1 line in the case of multiple cases matched to the same control or multiple controls matched to the same case. [Default: comparesheet.csv in /data]

Path to the genome reference
- Path to fasta reference [Default: ../data/Reference]

Path to bed exclude file
- Bed file containing regions to ignore [Default: Any .bed file in the /data directory]

Path to blacklisted fasta file
- Fasta sequences to ignore. [Default: Any *blacklist.gz file in the /data directory]

Type of structural variation
- The type of Structural Variants to detect: DELetions, INSertions, DUPlications, INVersions, BreakeND, or ALL types. [Default ALL]

Threads
- Number of threads to run on. Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples. [Default: max available]

### Auxilliary Parameters
- Minimum MAPQ quality for paired end mapping during SV discovery.  This is an integer value ranging from 0 to 255 and the optimal value can depend on the aligner used.  When in doubt, the user can plot the MAPQ distribution in a BAM file [Default: 1]
- Minimum distance on a reference for how far apart split-reads need to be [Default: 25 bp]
- Maximum distance on a reference for how far apart split-reads can be [Default: 40 bp]
- Minimum MAPQ quality for paired end mapping during SV genotyping. This is an integer value ranging from 0 to 255 and the optimal value can depend on the aligner used.  When in doubt, the user can plot the MAPQ distribution in a BAM file [Default: 5]
- gzipped output file for SV-reads
- Minimum SV site quality for merging [Default: 300]
- Minimum coverage for merging [Default: 10]
- Minimum SV length [Default: 0 bp]
- Maximum SV length [Default: 1,000,000 bp for merging; 500,000,000 bp for filtering]
- Only retain SVs with Precise tag during merging? This should only be used with whole exon data, not whole genome data [Default: no]
- Only retain SVs with Pass tag during merging? [Default: no for SV discovery; yes for SV genotyping]
- Minimum fraction of genotyped samples for filtering [Default: 0.75]
- Minimum mapping quality for CNVs [Default: 10]
- Ploidy [Default: 2]
- Minimum CNV size [Default: 1000]
- Window size for read-depth windows [Default: 10,000]
- Window offset for read-depth windows [Default: 10,000]
- Input BED file with the windows for read-depth windows
- Minimum fraction of the window that is callable for read-depth windows [Default: 0.25]
- Use mappable bases for window size? [Default: no]
- Scanning window size for GC fragment normalization [Default: 10,000]
- Scanning regions in BED format for GC fragment normalization
- Scan window selection for GC fragment normalization? [Default: no]
- Maximum CNV size [Default: 500,000,000]

## Outputs
- A folder for SV calling entitled ```/results/dellySV/``` containing
    - A subfolder ```/1calls``` containing SV calls for each case/control pair, as specified in the comparesheet, saved as a ```.bcf``` with an index
    - A subfolder ```/2prefilter``` containing filtered SV calls in a ```.bcf``` format with an index
    - A subfolder ```/3SV_filter``` containing a subfolder for each case sample. In each folder is a genotyped ```{sample_name}_geno.bcf``` file with index and the filtered genotypes in a ```{sample_name}_somatic.bcf``` file with index. 
- A folder for CNV calling entitled ```/results/dellyCNV/``` containing
    - A subfolder ```/1tumor``` containing CNV calls for each case sample, saved as a ```.bcf``` with an index. In addition, a ```_out.cov.gz``` file is saved for each sample to plot the read depth and copy number segmentation.  
    - A subfolder ```/2control``` containing a ```.bcf``` file with CNV calls for each control sample, with an index
    - A subfolder ```/3merge``` containing CNV calls for the matched case and control merged into a single ```.bcf``` file, with an index
    - A subfolder ```/4filter``` containing the CNV calls in the merged case control file that passed filtering in a ```.bcf``` format, with an index
    - A subfolder ```/5plot``` containing a bed file of the copy number segmentation, and a folder for each sample.  The folder contains a plot of the read depth of each chromosome and the entire genome for each sample. Copy number segmentation (after read number normalization) is laid on top of this to visualize the location of CNVs. 

## Cite

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.
DELLY: structural variant discovery by integrated paired-end and split-read analysis.
Bioinformatics. 2012 Sep 15;28(18):i333-i339.
https://doi.org/10.1093/bioinformatics/bts378

## Source

https://github.com/dellytools/delly

<hr>

[Code Ocean](https://codeocean.com/) is a cloud-based computational platform that aims to make it easy for researchers to share, discover, and run code.<br /><br />
[![Code Ocean Logo](images/CO_logo_68x36.png)](https://www.codeocean.com)