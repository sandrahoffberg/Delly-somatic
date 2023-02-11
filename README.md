## Germline structural variant & CNV calling in [Delly](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436805/)


Documentation: [Github](https://github.com/dellytools/delly)

Delly is an integrated structural variant (SV) prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read and long-read massively parallel sequencing data. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome.  This capsule does Structural Variant calling based on paired-end and split-reads and Copy Number Variation calling based on read depth in (Illumina) short reads only.

Delly includes many tunable parameters (e.g., for discovering variants with minimum quality, within specific window sizes, etc).  This capsule leaves most parameters as default.

Input: 
- BAM file that has been trimmed, aligned to reference, sorted, duplicate marked, read groups adjusted if necessary, indexed.  This capsule downloads data from an S3 bucket if a URL is provided in the App Panel. If no URL is provided, it will search for bam files in the ```/data``` folder.  It cannot use sample data from both locations in the same run. 
- Genome in fasta format (```.fa``` ending) with an index created in Samtools that ends in ```.fa.fai```  in the ```/data``` folder is required to identify split-reads
- Optional: mappability map, [downloaded here](https://gear.embl.de/data/delly/)
- Optional: - a BED file of repetitive regions to exclude (e.g., telomeres and centromeres) located in the ```/data``` folder.  The one used here was downloaded from [SpeedSeq Github page](https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed) and linked to from smoove Github page.
- Parameters provided in the App Panel
    - The URL of the S3 bucket with ```bam``` files, if they are not in the ```/data``` folder of this capsule. 
    - The location to downolad the ```bam``` files from S3 [Default: /results/data]
    - Whether to process samples individually or in batch for SV calling. Processing individually is better for high coverage samples, while processing in batch provides more power to detect variants for low coverage samples.  
    - Path to the genome reference
    - Path to the bed file of regions to exclude
    - Path to the genome annotation file 
    - Number of threads to run on. Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples. [Default: 1]


Output:
- Optional: the bam files, if downloaded from S3
- A folder for SV calling entitled ```/results/dellySV/``` containing
    - A subfolder ```/1calls``` containing SV calls for each sample, saved in a VCF with an index
    - A subfolder ```/2mergelist``` containing a bcf file with SV calls for each sample merged, called ```mergedSVdelly.bcf``` with an index
    - A subfolder ```/3geno``` containing more accurate SV genotypes for each sample, saved in a VCF with an index
    - A subfolder ```/4merge``` containing a bcf file with SV genotypes for each sample merged called ```merged.bcf``` with an index
    - If more than 20 samples present, a subfolder ```/5filter``` containing a filtered bcf file of the merged SV genotypes, called ```germline.bcf``` with an index
- A folder for CNV calling entitled ```/results/dellyCNV/``` containing
    - A subfolder ```/1CNVcall``` containing CNV calls for each sample, saved as a BCF with an index
    - A subfolder ```/2merge``` containing a vcf file with CNV calls for each sample merged, called ```dellyCNV.vcf``` with an index
    - A subfolder ```/3CNV``` containing more accurate CNV genotypes for each sample, saved in a BCF with an index.  In addition, a ```_out.cov.gz``` file is saved for each sample to plot the read depth and copy number segmentation. 
    - A subfolder ```/4merge``` containing a bcf file with SV genotypes for each sample merged called ```merged.bcf``` with an index
    - A subfolder ```/5filter``` containing a filtered bcf file of the merged SV genotypes, called ```filtered.bcf``` with an index
    - A subfolder ```/6plot``` containing 2 types of plots: 
        - If more than 100 samples present: a separate plot for each CNV is produced that show a histogram of how many samples have each number of copies.
        - Copy number segmentation.  A plot is made of the read depth of each chromosome and the entire genome for each sample. Copy number segmentation (after read number normalization) is laid on top of this to visualize the location of CNVs. 
