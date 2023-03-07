[![Code Ocean Logo](images/CO_logo_135x72.png)](http://codeocean.com/product)

<hr>

# Somatic structural variant & CNV calling in Delly

Delly is an integrated structural variant (SV) prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read and long-read massively parallel sequencing data. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome.  

This capsule does Structural Variant calling based on paired-end and split-reads and Copy Number Variation calling based on read depth in (Illumina) short reads from case (e.g., cancer) and matched control samples.

For more information, see the [Delly Github](https://github.com/dellytools/delly) page. 

## Input 
- BAM file that has been trimmed, aligned to reference, sorted, duplicate marked, read groups adjusted if necessary, indexed.  This capsule downloads data from an S3 bucket if a URL is provided in the App Panel. If no URL is provided, it will search for bam files in the **/data** directory.  It can use sample data from both locations in the same run. 
- In the **/data** directory, a genome in fasta format (```.fa``` ending) with an index created in Samtools that ends in ```.fa.fai```   is required to identify split-reads
- In the **/data** directory, a samplesheet, a tab-delimited (```.tsv```) sample description file where the first column is the sample ID (as in the VCF/BCF file) and the second column is either "tumor" or "control".
- In the **/data** directory, a comparesheet, a comma-delimited (```.csv```) file that has two samples on each line, the first is the control, the second is the matched tumor.  Sample names may appear on more than 1 line in the case of multiple cases matched to the same control or multiple controls matched to the same case. 
- Optional: In the **/data** directory, a mappability map, [downloaded here](https://gear.embl.de/data/delly/).  This is a ```fasta``` format file with repetitive regions of the genome hardmasked, and both an ```.fai``` and ```.gzi``` index. 
- Optional: In the **/data** directory, a BED file of repetitive regions to exclude (e.g., telomeres and centromeres).  The one used here was downloaded from [SpeedSeq Github page](https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed) and linked to from smoove Github page.



## Outputs
- Optional: the bam files, if downloaded from S3
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


## App Panel Parameters

### Main Parameters

- Path to samplesheet
- The URL of the S3 bucket with ```bam``` files, if they are not in the ```/data``` folder of this capsule. 
- The location to downolad the ```bam``` files from S3 [Default: /results/data]
- Path to comparesheet  
- Path to the genome reference
- Path to the bed file of regions to exclude
- Path to the blacklisted fasta file
- The type of Structural Variants to detect: DELetions, INSertions, DUPlications, INVersions, BreakeND, or ALL types. 
- Number of threads to run on. Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples. [Default: 1]

### Auxilliary Parameters

- Minimum paired end mapping quality [Default: 1]
- Minimum paired end quality for translocations [Default: 20]
- Minimum min. clipped length for a soft-clipped read [Default: 25]
- Minimum distance on a reference for how far apart split-reads need to be [Default: 25]
- Maximum distance on a reference for how far apart split-reads can be [Default: 40]
- Minimum mapping quality for genotyping [Default: 5]
- Optional: gzipped output file for SV-reads


## Citation

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.
DELLY: structural variant discovery by integrated paired-end and split-read analysis.
Bioinformatics. 2012 Sep 15;28(18):i333-i339.
https://doi.org/10.1093/bioinformatics/bts378


<hr>

[Code Ocean](https://codeocean.com/) is a cloud-based computational platform that aims to make it easy for researchers to share, discover, and run code.<br /><br />
[![Code Ocean Logo](images/CO_logo_68x36.png)](https://www.codeocean.com)