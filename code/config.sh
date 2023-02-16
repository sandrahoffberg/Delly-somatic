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



# URL of the s3 bucket with bams for analysis
if [ -z ${1} ]; then
    num_bam_files=$(find -L ../data -name "*.bam" | wc -l)
    if [ ${num_bam_files} -gt 0 ]; then
        echo "Using bams in the /data folder"
        data_dir="../data"
    else
        echo "There are no sample files available for this capsule.  Please either add some to the /data folder or provide a URL to an S3 bucket."
    fi
else
    s3_url=${1}
    
    # location to download the bams
    if [ -z ${2} ]; then
        data_dir=../results/data
    else
        data_dir=${2}
    fi
    
    #get BAM file(s) that have been trimmed, aligned to reference, sorted, read groups adjusted if necessary, and indexed
    echo "Downloading data from the S3 bucket at the provided URL."
    aws s3 sync --no-sign-request --only-show-errors ${s3_url} ${data_dir}
fi



if [ -z ${3} ]; then
  samplesheet=$(find -L ../data/ -name "samplesheet*.tsv")
else
  samplesheet=${3}
fi
# make the path absolute
samplesheet=$PWD/${samplesheet}


if [ -z ${4} ]; then
  comparesheet=$(find -L ../data/ -name "comparesheet*.csv")
else
  comparesheet=${4}
fi
# make the path absolute
comparesheet=$PWD/${comparesheet}





# If a reference file has not been specified in the app panel, find it
if [ -z ${5} ]; then
    reference=$(find -L ../data/ -name "*.fa")
    if [ "${reference}" = "" ]; then
        echo "Check your input data: there is no reference .fa file."
    fi
else
    reference=${5}
fi
# make the path absolute
reference=$PWD/${reference}


# If a bed file of regions to exclude has not been specified in the app panel, find it
if [ -z ${6} ]; then
  SVbed=$(find -L ../data/ -name "*.bed")
else
  SVbed=${6}
fi
# make the path absolute
SVbed=$PWD/${SVbed}


# If an annotation file has not been specified in the app panel, find it
if [ -z ${7} ]; then
    CNVmap=$(find -L ../data/ -name "*blacklist.gz")
else
    CNVmap=${7}
fi
# make the path absolute
CNVmap=$PWD/${CNVmap}



#Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples.
cores=$(get_cpu_count.py)

if [ -z ${8} ]; then
    num_threads=$((cores * 2))
else
    num_threads=${8}
fi

num_attached_bams=$(find -L ${data_dir} -name "*.bam" | wc -l)

min=$(( ${num_attached_bams} < ${num_threads} ? ${num_attached_bams} : ${num_threads} ))

export OMP_NUM_THREADS=${min}
