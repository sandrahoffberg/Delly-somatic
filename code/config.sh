#!/usr/bin/env bash

if [ $# -eq 0 ]; then
    echo "No arguments supplied"
else
    echo "args:"
    for i in $*; do 
        echo $i 
    done
    echo ""
fi


if [ -z ${1} ]; then
    samplesheet=$(find -L ../data/ -name "samplesheet*.tsv")
else
    samplesheet=${1}
fi

# Check that there is only 1 reference
number_samplesheet=$(find -L ../data/ -name "*samplesheet*" | wc -l )
if [ ${number_samplesheet} -eq 0 ]; then
    echo "Check your input data: there is no samplesheet."
    exit 1
elif [ ${number_samplesheet} -gt 1 ]; then
    echo "Check your input data: there are multiple samplesheets."
    exit 1
fi

# make the path absolute
samplesheet=$PWD/${samplesheet}


# URL of the s3 bucket with bams for analysis
if [ -z ${2} ]; then
    num_bam_files=$(find -L ../data -name "*.bam" | wc -l)
    if [ ${num_bam_files} -gt 0 ]; then
        echo "Using bams in the /data folder"
        data_dir="../data"
    else
        echo "There are no sample files available for this capsule.  Please either add some to the /data folder or provide a URL to an S3 bucket."
    fi
else
    s3_url=${2}
    
    # location to download the bams
    if [ -z ${3} ]; then
        data_dir=../scratch/data
    else
        data_dir=${3}
    fi
    
    #get BAM file(s) that have been trimmed, aligned to reference, sorted, read groups adjusted if necessary, and indexed
    echo "Downloading data from the S3 bucket at the provided URL."
    aws s3 sync --no-sign-request --only-show-errors ${s3_url} ${data_dir}
fi


if [ -z ${4} ]; then
    comparesheet=$(find -L ../data/ -name "comparesheet*.csv")
else
    comparesheet=${4}
fi

# Check that there is only 1 reference
number_comparesheet=$(find -L ../data/ -name "*compare*" | wc -l )
if [ ${number_comparesheet} -eq 0 ]; then
    echo "Check your input data: there is no comparesheet."
    exit 1
elif [ ${number_comparesheet} -gt 1 ]; then
    echo "Check your input data: there are multiple comparesheets."
    exit 1
fi

# make the path absolute
comparesheet=$PWD/${comparesheet}

# If a reference file has not been specified in the app panel, find it
if [ -z ${5} ]; then
    reference=$(find -L ../data/ -name "*.fa")
else
    reference=${5}
fi

# Check that there is only 1 reference
number_references=$(find -L ../data/ -name "*.fa" | wc -l )
if [ ${number_references} -eq 0 ]; then
    echo "Check your input data: there is no reference .fa file."
    exit 1
elif [ ${number_references} -gt 1 ]; then
    echo "Check your input data: there are multiple reference .fa files."
    exit 1
fi

# make the path absolute
reference=$PWD/${reference}


# If a bed file of regions to exclude has not been specified in the app panel, find it
if [ -z ${6} ]; then
  SVbed=$(find -L ../data/ -name "*.bed")
else
  SVbed=${6}
fi

# Check that there is only 1 bed file
number_SVbed=$(find -L ../data/ -name "*.bed" | wc -l )
if [ ${number_SVbed} -gt 1 ]; then
  echo "Check your input data: there are multiple bed files with regions to exclude."
  exit 1
fi

# make the path absolute
SVbed=$PWD/${SVbed}


# If an annotation file has not been specified in the app panel, find it
if [ -z ${7} ]; then
    CNVmap=$(find -L ../data/ -name "*blacklist.gz")
else
    CNVmap=${7}
fi

number_CNVmap=$(find -L ../data/ -name "*blacklist.gz" | wc -l )
if [ ${number_CNVmap} -gt 1 ]; then
    echo "Check your input data: there are multiple blacklist.gz files."
    exit 1
fi

# make the path absolute
CNVmap=$PWD/${CNVmap}


if [ -z ${8} ]; then
    svtype=ALL
else
    svtype="${8}"
fi


#Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples.

if [ -z ${9} ]; then
    num_threads=$(get_cpu_count.py)
else
    num_threads=${9}
fi

num_attached_bams=$(find -L ${data_dir} -name "*.bam" | wc -l)

# Check to make sure all of the necessary files are attached. 
if [ ${num_attached_bams} == 0 ]; then
    echo "Check your input data: there are no (indexed) bam files."
    exit 1
fi

min=$(( ${num_attached_bams} < ${num_threads} ? ${num_attached_bams} : ${num_threads} ))

export OMP_NUM_THREADS=${min}




## Auxilliary Parameters

## SV discovery options in delly call
if [ -z ${10} ]; then
    min_read_map_quality_discovery=
else
    min_read_map_quality_discovery=--map_qual ${10}
fi


if [ -z ${11} ]; then
    min_split_read_dist=
else
    min_split_read_dist=--minrefsep ${11}
fi


if [ -z ${12} ]; then
    max_split_read_dist=
else
    max_split_read_dist=--maxreadsep ${12}
fi


## SV genotyping options in delly call

if [ -z ${13} ]; then
    min_read_map_quality_genotyping=
else
    min_read_map_quality_genotyping=--geno_qual ${13}
fi


if [ -z ${14} ]; then
    SV_reads_output=
else
    SV_reads_output=--dump ${14}
fi


## merging options 

if [ -z ${15} ]; then
    min_SV_quality=
else
    min_SV_quality=--quality ${15}
fi


if [ -z ${16} ]; then
    min_coverage=
else
    min_coverage=--coverage ${16}
fi


if [ -z ${17} ]; then
    min_SV_size=
else
    min_SV_size=--min_size ${17}
fi


if [ -z ${18} ]; then
    max_SV_size=
else
    max_SV_size=--maxsize ${18}
fi


if [ ${19} == "no" ]; then
    precise_only=
else
    precise_only=--precise
fi


if [ ${20} == "no" ]; then
    pass_only=
else
    pass_only=--pass
fi



# delly filter options 

if [ -z ${21} ]; then
    min_genotyped_samples=
else
    min_genotyped_samples=--ratiogeno ${21}
fi


# delly CNV options

if [ -z ${22} ]; then
    min_mapping_quality_CNV=
else
    min_mapping_quality_CNV=--quality ${22}
fi


if [ -z ${23} ]; then
    ploidy=
else
    ploidy=--ploidy ${23}
fi


if [ -z ${24} ]; then
    min_CNV_size=
    min_CNV_size2=
else
    min_CNV_size=--cnv-size ${24}
    min_CNV_size2=--minsize ${24}
fi


# read depth window options
if [ -z ${25} ]; then
    window_size=
else
    window_size=--window-size ${25}
fi


if [ -z ${26} ]; then
    window_offset=
else
    window_offset=--window-offset ${26}
fi


if [ -z ${27} ]; then
    bed_intervals=
else
    bed_intervals=--bed-intervals ${27}
fi


if [ -z ${28} ]; then
    min_callable_window_fraction=
else
    min_callable_window_fraction=--fraction-window ${28}
fi


if [ ${29} == "no" ]; then
    mappable_bases_in_window=
else
    mappable_bases_in_window=--adaptive-windowing
fi


# GC fragment normalization options
if [ -z ${30} ]; then
    scan_window_size=
else
    scan_window_size=--scan-windowing ${30}
fi


if [ -z ${31} ]; then
    scanning_regions_in_BED=
else
    scanning_regions_in_BED=--scan-regions ${31}
fi


if [ ${32} == "no" ]; then
    scan_window=
else
    scan_window=--no-window-selection
fi


# Delly classify options

if [ -z ${33} ]; then
    max_CNV_size=
else
    max_CNV_size=--maxsize ${33}
fi
